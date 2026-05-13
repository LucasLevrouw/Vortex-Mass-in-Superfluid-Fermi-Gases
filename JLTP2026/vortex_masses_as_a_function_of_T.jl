### Code written by Lucas Levrouw (2026)
### This script generates figure 4 of https://arxiv.org/abs/2512.22099
### If using this code, please cite the above paper.



# Run calculate-profiles.jl first to generate the data

using Revise
using VortexMass; const PROJECT_ROOT = pkgdir(VortexMass)
using DelimitedFiles
using NumericalIntegration

# Create necessary directories
mkpath("$PROJECT_ROOT/JLTP2026/tables/")
mkpath("$PROJECT_ROOT/JLTP2026/figures/")

datapath = "$PROJECT_ROOT/data/"

read_from_table = false  # set to true to read masses from saved tables

zeta_crit_data = readdlm(normpath(joinpath(PROJECT_ROOT, "JLTP2026/critical_ζ.dat")), '\t')
zeta_crit_interactions = zeta_crit_data[:,1]
zeta_crit_vals = zeta_crit_data[:,2]


temperatures = 0.0:0.05:0.95
interactions = -1.0:1.0:1.0
rel_zetaVals = [0.0, 0.5, 0.7, 0.8, 0.9]
zetaCrit = zeros(length(interactions))

bc = "neumann"

L_kF = 150

# Pre-compute zetaCrit from the lookup table (doesn't require profiles)
for j in eachindex(interactions)
    idx = findfirst(interactions[j] .== zeta_crit_interactions)
    if idx !== nothing
        zetaCrit[j] = zeta_crit_vals[idx]
    end
end

if !read_from_table
    println("Loading profile data...")
    r_array  = Array{Union{Vector{Union{Float64, Missing}}, Missing}}(missing, length(temperatures), length(interactions), length(rel_zetaVals))
    relaxed_profiles  = copy(r_array)

    for k in eachindex(rel_zetaVals), j in eachindex(interactions), i in eachindex(temperatures)
        int = interactions[j]
        T = temperatures[i]
        β = (T == 0) ? 1000 : 1/(T*find_Tc(int))

        ζ = rel_zetaVals[k] * zetaCrit[j]

        filepath = normpath("$datapath/relaxed_profile_$(bc)_int=$(int)_T=$(T)Tc_ζ=$(rel_zetaVals[k])ζc.csv")

        if !isfile(filepath)
            Δ, _ = find_Δ_μ(β, int, ζ)
            if Δ ≠ 0.0
                println("Profile not found with nonzero Δ = $Δ at int = $int, T/T_c = $T, ζ/ζ_c = $(rel_zetaVals[k]), skipping...")
            end
            continue
        end

        relaxed_data = readdlm(filepath,',')
        r_array[i,j,k] = relaxed_data[:,1]
        relaxed_profiles[i,j,k] = relaxed_data[:,2]
    end

    associated_masses = Array{Union{Float64, Missing}}(missing, length(temperatures), length(interactions), length(rel_zetaVals))
    internal_masses = copy(associated_masses)
    imbalance_masses = copy(associated_masses)

    associated_masses_tail = copy(associated_masses)
    internal_masses_tail = copy(associated_masses)
    imbalance_masses_tail = copy(associated_masses)

    mass_scales = copy(associated_masses)
    healing_lengths = copy(associated_masses)
    coefficients_internal = copy(associated_masses)
    correction_factor_a = copy(associated_masses)
    correction_factor_i = copy(associated_masses)

    println("Calculating masses...")

    for k in eachindex(rel_zetaVals), j in eachindex(interactions), i in eachindex(temperatures)
        int = interactions[j]
        T = temperatures[i]
        Tc = find_Tc(int)
        β = (T == 0) ? 1000 : 1/(T*Tc)
        ζ = zetaCrit[j] * rel_zetaVals[k]
        (; Δ, μ, A_itp, Dt_itp, Ct, Et, Q, Rt, G, ξ2) = calculate_mf_eft_parameters(β, int, ζ, warn = false)
        if Δ == 0.0
            continue
        end
        if ismissing(relaxed_profiles[i,j,k]) || ismissing(r_array[i,j,k])
            println("Missing profile for int = $int, T = $T Tc, ζ = $ζ ζc")
            continue
        end
        u = Δ * relaxed_profiles[i,j,k]

        healing_lengths[i,j,k] = ξ2
        local L = L_kF / ξ2

        r = r_array[i,j,k]

        ρ_tot = 1/2 * broadcast(Φ -> numberDensity_mf(β, μ, ζ, Φ)[1], u)
        ρ_s = @. 2 * Ct * u^2
        ρ_n = ρ_tot - ρ_s
        ρ_imb = 1/2 * broadcast(Φ -> imbalanceDensity_mf(β, μ, ζ, Φ)[1], u)

        ρ_∞ = 1/6π^2
        ρ_s_∞ = 2 * Ct * Δ^2
        ρ_n_∞ = ρ_∞ - ρ_s_∞
        ρ_imb_∞ = 1/2 * imbalanceDensity_mf(β, μ, ζ, Δ)[1]

        δρ_n_∞ = -Δ^2 * Dt_int(β, μ, ζ, Δ)[1] + ρ_s_∞

        coefficients_internal[i,j,k] = δρ_n_∞

        idx = findfirst(x -> x ≥ L * ξ2, r)

        m_a = 2π * integrate(r[1:idx], (r .* (ρ_s_∞ .- ρ_s))[1:idx])
        m_i = 2π * integrate(r[1:idx], (r .* (ρ_n .- ρ_n_∞))[1:idx])
        m_imb = 2π * integrate(r[1:idx], (r .* (ρ_imb .- ρ_imb_∞))[1:idx])

        associated_masses[i,j,k] = m_a
        internal_masses[i,j,k] = m_i
        imbalance_masses[i,j,k] = m_imb

        cutoff_radius = 5
        idx_cutoff = findfirst(x -> x ≥ cutoff_radius * ξ2, r)

        m_a_tail = idx_cutoff > idx ? missing : 2π * integrate(r[idx_cutoff:idx], (r .* (ρ_s_∞ .- ρ_s))[idx_cutoff:idx])
        m_i_tail = idx_cutoff > idx ? missing : 2π * integrate(r[idx_cutoff:idx], (r .* (ρ_n .- ρ_n_∞))[idx_cutoff:idx])
        m_imb_tail = idx_cutoff > idx ? missing : 2π * integrate(r[idx_cutoff:idx], (r .* (ρ_imb .- ρ_imb_∞))[idx_cutoff:idx])

        associated_masses_tail[i,j,k] = m_a_tail
        internal_masses_tail[i,j,k] = m_i_tail
        imbalance_masses_tail[i,j,k] = m_imb_tail

        mass_scales[i,j,k] = π * ξ2^2 * ρ_s_∞

        m_a0 = 2π * integrate(r, (r .* (ρ_s_∞ .- ρ_s))) + π * ξ2^2 * ρ_s_∞ * log(ξ2/r[end])
        m_i0 = 2π * integrate(r, (r .* (ρ_n .- ρ_n_∞))) + π * ξ2^2 * δρ_n_∞ * log(ξ2/r[end])

        correction_factor_a[i,j,k] = exp(-m_a0 / (π * ξ2^2 * ρ_s_∞))
        correction_factor_i[i,j,k] = exp(-m_i0 / (π * ξ2^2 * δρ_n_∞))
    end

    total_masses = associated_masses .+ internal_masses

    open("$PROJECT_ROOT/JLTP2026/tables/vortex_masses_as_function_of_T.csv", "w") do io
        write(io, "# Vortex masses as a function of temperature, for selected interactions and imbalance values\n")
        write(io, "# Generated by JLTP2026/calculate-and-plot-vortex-masses-as-a-function-of-T.jl\n")
        write(io, "# Columns: T/T_c, (k_F a_s)^(-1), zeta/zeta_c, M_a/(m k_F), M_i/(m k_F), M_imb/(m k_F), alpha_a, alpha_i\n")
        for k in eachindex(rel_zetaVals), j in eachindex(interactions), i in eachindex(temperatures)
            writedlm(io, [temperatures[i] interactions[j] rel_zetaVals[k] 2*coalesce(associated_masses[i,j,k], NaN) 2*coalesce(internal_masses[i,j,k], NaN) 2*coalesce(imbalance_masses[i,j,k], NaN) coalesce(correction_factor_a[i,j,k], NaN) coalesce(correction_factor_i[i,j,k], NaN)], ',')
        end
    end
else
    data = readdlm("$PROJECT_ROOT/JLTP2026/tables/vortex_masses_as_function_of_T.csv", ',', comments=true)
    nT, ni, nz = length(temperatures), length(interactions), length(rel_zetaVals)
    associated_masses   = reshape(data[:, 4], nT, ni, nz) ./ 2
    internal_masses     = reshape(data[:, 5], nT, ni, nz) ./ 2
    imbalance_masses    = reshape(data[:, 6], nT, ni, nz) ./ 2
    correction_factor_a = reshape(data[:, 7], nT, ni, nz)
    correction_factor_i = reshape(data[:, 8], nT, ni, nz)
    total_masses        = associated_masses .+ internal_masses
end

# Plots

include("initialize_makie.jl")


zeta_colors = [ColorSchemes.tab10[1],ColorSchemes.tab10[4],ColorSchemes.tab10[2],ColorSchemes.tab10[5],ColorSchemes.tab10[3]]
zeta_markers = [:rect, :diamond, :utriangle, :dtriangle, :cross]

R_colors = [ColorSchemes.tab10[7],ColorSchemes.tab10[9],ColorSchemes.tab10[8]]



### Vortex mass in function of T
fig = Figure(size = (490pt,196pt), figure_padding = (0,7,0.5pt,0.5pt))
axes = [Axis(fig[1, j], xlabel = L"T/T_c", ylabel = L"\frac{M_{\text{tot}}}{M_{\text{tot}}(\zeta = T = 0)}", xticks = 0.0:0.2:1.0, limits = ((-0.05,1.0),(nothing,nothing))) for j in eachindex(interactions)]
linkaxes!(axes...)

hideydecorations!.(axes[2], grid=false)
hideydecorations!.(axes[3], grid=false)

markers = [:+, :circle, :rect, :utriangle, :star5, :diamond, :dtriangle]

missing_to_nan(x) = ismissing(x) ? NaN : x


for j in eachindex(interactions), k in eachindex(rel_zetaVals)
    if all(ismissing.(total_masses[:,j,k]))
        continue
    end
    scatter!(axes[j], temperatures, missing_to_nan.(total_masses[:,j, k])./total_masses[begin, j, begin], marker = zeta_markers[k],
    color = zeta_colors[k],
    label=L"ζ/ζ_c=%$(rel_zetaVals[k])")
end

for j in eachindex(interactions)
    Label(fig[1,j,Top()], L"(k_F a_s)^{-1} = %$(Int64(interactions[j]))", halign = :center, padding = (0,0,4,0))
end

Legend(fig[0,1:3], axes[1], orientation = :horizontal)

save("$(PROJECT_ROOT)/JLTP2026/figures/fig4_vortex_mass_as_a_function_of_T.pdf", fig)
