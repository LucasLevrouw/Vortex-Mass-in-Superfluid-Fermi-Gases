### Code written by Lucas Levrouw (2026)
### This script generates figure 7 of https://arxiv.org/abs/2512.22099
### If using this code, please cite the above paper.

# Run calculate-profiles.jl first to generate the data


using Revise
using VortexMass; const PROJECT_ROOT = pkgdir(VortexMass)
using OrdinaryDiffEq
using DelimitedFiles
using Trapz

# Create necessary directories
mkpath("$PROJECT_ROOT/JLTP2026/tables/")
mkpath("$PROJECT_ROOT/JLTP2026/figures/")

read_from_table = false  # set to true to read from saved tables instead of profiles

zeta_crit_data = readdlm(normpath(joinpath(PROJECT_ROOT, "JLTP2026/critical_ζ.dat")), '\t')
zeta_crit_interactions = zeta_crit_data[:,1]
zeta_crit_vals = zeta_crit_data[:,2]


temperatures = [0.0, 0.2, 0.4, 0.6]
interactions = -1.0:1.0:1.0
rel_zetaVals = [0.8,0.9]
zetaCrit = zeros(length(interactions))

bc = "neumann"

# Pre-compute zetaCrit from the lookup table
for j in eachindex(interactions)
    idx = findfirst(interactions[j] .== zeta_crit_interactions)
    if idx !== nothing
        zetaCrit[j] = zeta_crit_vals[idx]
    end
end

# Storage arrays for profile-derived quantities
r_xi_data        = Array{Union{Vector{Float64}, Missing}}(missing, length(temperatures), length(interactions), length(rel_zetaVals))
rho_tot_norm_data = Array{Union{Vector{Float64}, Missing}}(missing, length(temperatures), length(interactions), length(rel_zetaVals))
rho_s_norm_data  = Array{Union{Vector{Float64}, Missing}}(missing, length(temperatures), length(interactions), length(rel_zetaVals))
rho_n_norm_data  = Array{Union{Vector{Float64}, Missing}}(missing, length(temperatures), length(interactions), length(rel_zetaVals))
rho_imb_norm_data = Array{Union{Vector{Float64}, Missing}}(missing, length(temperatures), length(interactions), length(rel_zetaVals))

if !read_from_table
    for k in eachindex(rel_zetaVals), j in eachindex(interactions), i in eachindex(temperatures)
        int = interactions[j]
        T = temperatures[i]
        β = (T == 0) ? 1000 : 1/(T*find_Tc(int))
        ζ = rel_zetaVals[k] * zetaCrit[j]
        println("int = $int, T/T_F = $(1/β), ζ = $ζ")

        (; Δ, μ, A_itp, Dt_itp, Ct, Et, Q, Rt, G, ξ2) = calculate_mf_eft_parameters(β,int,ζ)

        filepath = normpath(joinpath(PROJECT_ROOT,"data/relaxed_profile_$(bc)_int=$(int)_T=$(T)Tc_ζ=$(rel_zetaVals[k])ζc.csv"))
        if !isfile(filepath)
            println("No profile data for int = $int, T/T_F = $(1/β), ζ = $ζ, skipping...")
            continue
        end
        relaxed_data = readdlm(filepath,',')

        u = Δ .* relaxed_data[:,2]
        r = relaxed_data[:,1]
        ξ2 = sqrt(2*Ct /(G*Δ^2))

        ρ_∞ = 1/6π^2
        ρ_tot = 1/2*broadcast(Φ->numberDensity_mf(β, μ, ζ, Φ)[1], u)
        ρ_s   = @. 2*Ct*u^2
        ρ_n   = ρ_tot - ρ_s
        Δρ    = 1/2*broadcast(Φ->imbalanceDensity_mf(β, μ, ζ, Φ)[1], u)

        r_xi_data[i,j,k]        = r/ξ2
        rho_tot_norm_data[i,j,k] = ρ_tot/ρ_∞
        rho_s_norm_data[i,j,k]  = ρ_s/ρ_∞
        rho_n_norm_data[i,j,k]  = ρ_n/ρ_∞
        rho_imb_norm_data[i,j,k] = Δρ/ρ_∞

        open("$PROJECT_ROOT/JLTP2026/tables/fig7_densities_finT_T=$(T)_int=$(int)_zeta=$(rel_zetaVals[k])zetac.csv", "w") do io
            write(io, "# This data corresponds to (a subplot of) figure 7 of https://arxiv.org/abs/2512.22099\n# If using this data, please cite the above paper.\n")
            writedlm(io, permutedims(["r/xi", "rho_tot/rho_inf", "rho_s/rho_inf", "rho_n/rho_inf", "rho_imbalance/rho_inf"]), ',')
            writedlm(io, hcat(r/ξ2, ρ_tot/ρ_∞, ρ_s/ρ_∞, ρ_n/ρ_∞, Δρ/ρ_∞), ',')
        end
    end
else
    for k in eachindex(rel_zetaVals), j in eachindex(interactions), i in eachindex(temperatures)
        int = interactions[j]
        T = temperatures[i]
        tbl = readdlm("$PROJECT_ROOT/JLTP2026/tables/fig7_densities_finT_T=$(T)_int=$(int)_zeta=$(rel_zetaVals[k])zetac.csv", ',', comments=true, skipstart=3)
        r_xi_data[i,j,k]        = Float64.(tbl[:,1])
        rho_tot_norm_data[i,j,k] = Float64.(tbl[:,2])
        rho_s_norm_data[i,j,k]  = Float64.(tbl[:,3])
        rho_n_norm_data[i,j,k]  = Float64.(tbl[:,4])
        rho_imb_norm_data[i,j,k] = Float64.(tbl[:,5])
    end
end

include("initialize_makie.jl")

k_fixed = 2
rmax = 3

rel_zeta = rel_zetaVals[k_fixed]
println(L"ζ/ζ_c = %$(rel_zeta)")


fig = Figure(size = (490pt,350pt), figure_padding = (0,1pt,1pt,0))

ax = [Axis(fig[i,j], xlabel=L"r/\xi", ylabel=L"\rho/\rho_{\text{tot},\infty}") for i in eachindex(temperatures), j in eachindex(interactions)]

linewidths = [0.75pt, 1.125pt, 1.5pt, 0.5pt]

for i in eachindex(temperatures), j in eachindex(interactions)
    k = k_fixed
    if i ≠ length(temperatures)
        hidexdecorations!(ax[i,j], grid=false)
    end
    if j ≠ 1
        hideydecorations!(ax[i,j], grid=false)
    end

    int = interactions[j]
    ζ = rel_zetaVals[k] * zetaCrit[j]
    T = temperatures[i]*find_Tc(int)
    β = (T == 0) ? 1e3 : 1/T
    Δ, μ = find_Δ_μ(β, int, ζ)

    if Δ ≠ 0.0 && !ismissing(r_xi_data[i,j,k])
        (; Δ, μ, A_itp, Dt_itp, Ct, Et, Q, Rt, G, ξ2) = calculate_mf_eft_parameters(β,int,ζ)

        ρ_∞ = 1/6π^2
        ρ_s∞ = 2*Ct*Δ^2
        ρ_n∞ = ρ_∞ - ρ_s∞
        δρ_n∞ = ρ_s∞ - Δ^2 * Dt_int(β, μ, ζ, Δ)[1]

        r_xi        = r_xi_data[i,j,k]
        rho_tot_norm = rho_tot_norm_data[i,j,k]
        rho_s_norm  = rho_s_norm_data[i,j,k]
        rho_n_norm  = rho_n_norm_data[i,j,k]
        rho_imb_norm = rho_imb_norm_data[i,j,k]

        xlims!(ax[i,j], 0, rmax)
        ylims!(ax[i,j], -0.05, 1.05)

        lines!(ax[i,j], r_xi, rho_tot_norm, label="total", color=mass_colors[1], linewidth=linewidths[1])
        lines!(ax[i,j], r_xi, rho_s_norm, label="superfluid", color=mass_colors[2], linewidth=linewidths[2])
        lines!(ax[i,j], r_xi, rho_n_norm, label="normal", color=mass_colors[3], linewidth=linewidths[3])
        lines!(ax[i,j], r_xi, rho_imb_norm, label="imbalance", color = :black, linestyle = :dot)

        nanif(x,bool) = bool ? NaN : x

        lines!(ax[i,j], r_xi, nanif.(ρ_s∞/ ρ_∞ .* (1 .- 1/2 * r_xi.^-2), r_xi .<= 1), color = asymptotic_colors[2], linestyle=:dash, alpha = 0.8)
        lines!(ax[i,j], r_xi, nanif.(ρ_n∞ / ρ_∞ .+ δρ_n∞ / ρ_∞ * 1/2 * r_xi.^-2, r_xi .<= 1), color = asymptotic_colors[3], linestyle=:dash , alpha =0.8)
    end
end

elem_t = LineElement(color = mass_colors[1], linewidth = linewidths[1])
elem_s = [LineElement(color = mass_colors[2], linewidth = linewidths[2], points = Point2f[(0, 0.6), (1, 0.6)]), LineElement(color = asymptotic_colors[2], points = Point2f[(0, 0.4), (1, 0.4)], linestyle=:dash, alpha = 0.8)]
elem_n = [LineElement(color = mass_colors[3], linewidth = linewidths[3], points = Point2f[(0, 0.6), (1, 0.6)]), LineElement(color = asymptotic_colors[3], points = Point2f[(0, 0.4), (1, 0.4)], linestyle=:dash, alpha = 0.8)]
elem_i = LineElement(color = :black, linewidth = linewidths[3], linestyle = :dot)
Legend(fig[0, :], [elem_t, elem_s, elem_n, elem_i], ["total", "superfluid", "normal", "imbalance"], orientation = :horizontal)


rowgap!(fig.layout, 10)

heights = [0.15, 0.33, 0.06, 0.05]

for i in eachindex(temperatures), j in eachindex(interactions)
    label = L"T/T_c = %$(temperatures[i])"
    text!(
        ax[i,j], 0.98, j==1 && i==3 ? 0.3 : heights[i],
        text = label,
        align = (:right, :bottom),
        space = :relative,
        fontsize = 9pt,
    )
end

for j in eachindex(interactions)
    Label(fig[1,j,Top()],  L"(k_F a_s)^{-1} = %$(interactions[j])", halign = :center, padding = (0,0,4,0))
end

Box(fig[0,1], color = :white, width = 45pt, tellwidth = false, halign = :left)
Label(fig[0,1], L"\zeta = %$(rel_zetaVals[k_fixed])\zeta_c", halign = :left, valign = :center, tellwidth = false, fontsize = 9pt, padding = (6pt,0,0,0))

save("$PROJECT_ROOT/JLTP2026/figures/fig7_densities_finT_ζ=$(rel_zetaVals[k_fixed])ζc.pdf", fig)

fig
