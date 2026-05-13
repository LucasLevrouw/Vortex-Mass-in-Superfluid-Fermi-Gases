using Test
using Random
using VortexMass
using QuadGK
using NonlinearSolve
using SciMLBase
using BSplineKit

import VortexMass: xi, En, En_minus_xi, X, f1, f2, f3, f4, I5, I6

include(joinpath(@__DIR__, "eft_coefficients_and_mean_field_old.jl"))

# Fixed parameters
const β_test = 3.0
const μ_test = 0.5
const ζ_test = 0.2
const Δ_test = 0.4
const int_test = 0.5

@testset "EFT coefficients: current vs. reference implementations" begin

    @testset "Helper function X (ϵ > ζ)" begin
        # At ϵ > ζ, X_old and X use identical formulas; results must agree to machine precision.
        ϵ = 0.8  # > ζ_test = 0.2
        @test X_old(ϵ, β_test, ζ_test) ≈ X(ϵ, β_test, ζ_test)  rtol=1e-12
    end

    @testset "Helper functions f1–f4 (ϵ > ζ)" begin
        # f1_old uses X_old; at ϵ > ζ both X_old and X evaluate the same expression.
        # f2_old through f4_old use inline formulas identical to the current implementations.
        ϵ = 0.8
        @test f1_old(ϵ, β_test, ζ_test) ≈ f1(ϵ, β_test, ζ_test)  rtol=1e-12
        @test f2_old(ϵ, β_test, ζ_test) ≈ f2(ϵ, β_test, ζ_test)  rtol=1e-12
        @test f3_old(ϵ, β_test, ζ_test) ≈ f3(ϵ, β_test, ζ_test)  rtol=1e-12
        @test f4_old(ϵ, β_test, ζ_test) ≈ f4(ϵ, β_test, ζ_test)  rtol=1e-12
    end

    @testset "EFT action coefficients" begin
        # C_int_old, Dt_int_old, etc. use the current f2/f3/f4 from VortexMass internally,
        # so agreement should be exact up to quadrature tolerance.
        @test C_int_old(β_test, μ_test, ζ_test, Δ_test)[1]  ≈ C_int(β_test, μ_test, ζ_test, Δ_test)[1]   rtol=1e-8
        @test Dt_int_old(β_test, μ_test, ζ_test, Δ_test)[1] ≈ Dt_int(β_test, μ_test, ζ_test, Δ_test)[1]  rtol=1e-8
        @test Et_int_old(β_test, μ_test, ζ_test, Δ_test)[1] ≈ Et_int(β_test, μ_test, ζ_test, Δ_test)[1]  rtol=1e-8
        @test Q_int_old(β_test, μ_test, ζ_test, Δ_test)[1]  ≈ Q_int(β_test, μ_test, ζ_test, Δ_test)[1]   rtol=1e-6
        @test Rt_int_old(β_test, μ_test, ζ_test, Δ_test)[1] ≈ Rt_int(β_test, μ_test, ζ_test, Δ_test)[1]  rtol=1e-6
        @test G_int_old(β_test, μ_test, ζ_test, Δ_test)[1]  ≈ G_int(β_test, μ_test, ζ_test, Δ_test)[1]   rtol=1e-8
    end

    @testset "Gap equation and number density" begin
        # A_int_old and numberDensity_mf_old use exp(-β*En)*cosh(β*ζ) without case splitting;
        # for the chosen parameters (En > ζ throughout the dominant integration region) the
        # values should agree with the current numerically-stabilised implementations.
        @test A_int_old(β_test, μ_test, ζ_test, Δ_test, int_test)[1]       ≈ A_int(β_test, μ_test, ζ_test, Δ_test, int_test)[1]       rtol=1e-8
        @test numberDensity_mf_old(β_test, μ_test, ζ_test, Δ_test)[1]      ≈ numberDensity_mf(β_test, μ_test, ζ_test, Δ_test)[1]      rtol=1e-8
        @test Omega_sp_old(β_test, μ_test, ζ_test, Δ_test, int_test)[1]    ≈ Omega_sp(β_test, μ_test, ζ_test, Δ_test, int_test)[1]    rtol=1e-6
    end

    @testset "Mean-field saddle point" begin
        # find_Δ_μ_old returns a raw NonlinearSolve solution; sol[1] = Δ, sol[2] = μ.
        sol_old = find_Δ_μ_old(β_test, int_test)
        Δ_new, μ_new = find_Δ_μ(β_test, int_test)
        @test sol_old[1] ≈ Δ_new  rtol=1e-6
        @test sol_old[2] ≈ μ_new  rtol=1e-6
    end
end

# Randomized tests: sample N_trials parameter sets from physically reasonable ranges and check
# that old and new implementations agree.  The seed is fixed so results are reproducible.
#
# Parameter constraints:
#   β*ζ ≤ 5*0.4 = 2.0 — keeps exp(-β*ϵ)*cosh(β*ζ) well within Float64 range for X_old.
#   Δ > 0             — gap must be positive for En and the coefficient integrals to be well-defined.
#   ζ < Δ             — physical FFLO stability condition; ensures the saddle-point integrals
#                        (which have En(k,μ,Δ) ≥ Δ > ζ) avoid the ϵ < ζ regime entirely in A_int/numberDensity.
#
# Helper functions (X, f1–f4) are tested at both ϵ > ζ and ϵ < ζ.  In the ϵ < ζ branch the
# old and new formulations are mathematically equal but take different floating-point paths;
# the tolerance is relaxed slightly to 1e-10 to allow for rounding differences.

const N_TRIALS = 5

@testset "EFT coefficients: randomized parameters (N=$N_TRIALS)" begin
    seed = rand(Random.RandomDevice(), UInt)
    @info "Randomized EFT coefficient tests using seed=$seed (use Random.MersenneTwister($seed) to reproduce)"
    rng = Random.MersenneTwister(seed)

    for trial in 1:N_TRIALS
        β   = 1.0 + 4.0  * rand(rng)   # ∈ [1, 5]
        μ   = -0.5 + 1.5 * rand(rng)   # ∈ [-0.5, 1.0]
        ζ   = 0.05 + 0.35 * rand(rng)  # ∈ [0.05, 0.4]  →  β*ζ ≤ 2.0
        Δ   = ζ + 0.1 + 0.5 * rand(rng) # > ζ, ∈ [ζ+0.1, ζ+0.6]
        int = -1.0 + 3.0 * rand(rng)   # ∈ [-1, 2]

        ϵ_hi = ζ + 0.05 + 0.5 * rand(rng)   # > ζ
        ϵ_lo = ζ * rand(rng) * 0.9           # ∈ (0, ζ), strictly less than ζ

        @testset "trial $trial (β=$(round(β,digits=2)), μ=$(round(μ,digits=2)), ζ=$(round(ζ,digits=2)), Δ=$(round(Δ,digits=2)))" begin

            @testset "X and f1–f4, ϵ > ζ" begin
                @test X_old(ϵ_hi, β, ζ)  ≈ X(ϵ_hi, β, ζ)   rtol=1e-10
                @test f1_old(ϵ_hi, β, ζ) ≈ f1(ϵ_hi, β, ζ)  rtol=1e-10
                @test f2_old(ϵ_hi, β, ζ) ≈ f2(ϵ_hi, β, ζ)  rtol=1e-10
                @test f3_old(ϵ_hi, β, ζ) ≈ f3(ϵ_hi, β, ζ)  rtol=1e-10
                @test f4_old(ϵ_hi, β, ζ) ≈ f4(ϵ_hi, β, ζ)  rtol=1e-10 # sometime errors! Check!
            end

            @testset "X and f1, ϵ < ζ" begin
                # Different floating-point paths; mathematically equivalent for small β*ζ.
                @test X_old(ϵ_lo, β, ζ)  ≈ X(ϵ_lo, β, ζ)   rtol=1e-10
                @test f1_old(ϵ_lo, β, ζ) ≈ f1(ϵ_lo, β, ζ)  rtol=1e-10
            end

            @testset "EFT action coefficients" begin
                @test C_int_old(β, μ, ζ, Δ)[1]  ≈ C_int(β, μ, ζ, Δ)[1]   rtol=1e-7
                @test Dt_int_old(β, μ, ζ, Δ)[1] ≈ Dt_int(β, μ, ζ, Δ)[1]  rtol=1e-7
                @test Et_int_old(β, μ, ζ, Δ)[1] ≈ Et_int(β, μ, ζ, Δ)[1]  rtol=1e-7
                @test Q_int_old(β, μ, ζ, Δ)[1]  ≈ Q_int(β, μ, ζ, Δ)[1]   rtol=1e-5
                @test Rt_int_old(β, μ, ζ, Δ)[1] ≈ Rt_int(β, μ, ζ, Δ)[1]  rtol=1e-5
                @test G_int_old(β, μ, ζ, Δ)[1]  ≈ G_int(β, μ, ζ, Δ)[1]   rtol=1e-7
            end

            @testset "Gap equation and number density" begin
                @test A_int_old(β, μ, ζ, Δ, int)[1]       ≈ A_int(β, μ, ζ, Δ, int)[1]       rtol=1e-7
                @test numberDensity_mf_old(β, μ, ζ, Δ)[1] ≈ numberDensity_mf(β, μ, ζ, Δ)[1] rtol=1e-7
                @test Omega_sp_old(β, μ, ζ, Δ, int)[1]    ≈ Omega_sp(β, μ, ζ, Δ, int)[1]    rtol=1e-5
            end
        end
    end
end
