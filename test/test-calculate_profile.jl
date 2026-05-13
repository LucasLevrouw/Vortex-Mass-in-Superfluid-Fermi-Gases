using Revise
using VortexMass; const PROJECT_ROOT = pkgdir(VortexMass)
using OrdinaryDiffEq
using Optim
using DelimitedFiles


T = 0.0
int = 3.0
ζ = 0.0

R_ξ2 = 100
N = 1000
l = 1

Tc = find_Tc(int)

β = (T == 0) ? 1000 : 1/(T*Tc)

println("Starting test for (k_F a_s)^-1 = $int, T/T_F = $(1/β)")

(; ξ2) = calculate_mf_eft_parameters(β, int, ζ)

R = R_ξ2 * ξ2

t = 1e4

output_dir = joinpath(PROJECT_ROOT, "test", "output")
mkpath(output_dir)

bc_variants = [
    (:neumann,           "neumann"),
    (:dirichlet,         "dirichlet"),
    ((:dirichlet, 0.95), "dirichlet_0.95"),
    ((:neumann, 0.0),    "neumann_0.0"),
]

diffs = Dict{String, Float64}()

for (bc, label) in bc_variants
    println("\nRunning bc = $bc ...")
    @time profile = calculate_profile(β, int, ζ, R, N, t, bc)

    outfile = joinpath(output_dir, "profile_$(label).csv")

    if isfile(outfile)
        profile_saved = readdlm(outfile)
        max_diff = maximum(abs, profile - profile_saved)
        diffs[label] = max_diff
        status = max_diff < 1e-4 ? "PASS" : "FAIL"
        println("  $status  bc=$bc  (max diff = $max_diff)")
    else
        println("  Saving baseline for bc=$bc → $outfile")
    end

    writedlm(outfile, profile)
end

if isempty(diffs)
    println("\nBaselines saved — rerun to compare.")
elseif all(v < 1e-4 for v in values(diffs))
    println("\nAll profiles match saved baselines (tol = 1e-4).")
else
    println("\nSome profiles do NOT match — check the diffs above.")
end
