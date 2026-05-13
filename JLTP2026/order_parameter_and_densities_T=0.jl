### Code written by Lucas Levrouw (2026)
### This script generates figures 1 and 6 of https://arxiv.org/abs/2512.22099
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


temperatures = 0.0 #[0.0, 0.5, 0.9]
interactions = -1.0:1.0:1.0
rel_zetaVals = [0.0,0.7,0.8, 0.9]
zetaCrit = zeros(length(interactions))

bc = "neumann"

# Pre-compute zetaCrit from the lookup table
for j in eachindex(interactions)
    idx = findfirst(interactions[j] .== zeta_crit_interactions)
    if idx !== nothing
        zetaCrit[j] = zeta_crit_vals[idx]
    end
end

# Storage arrays for fig 6 (order parameter and current)
r_xi_fig6      = Array{Union{Vector{Float64}, Missing}}(missing, length(interactions), length(rel_zetaVals))
f_data         = Array{Union{Vector{Float64}, Missing}}(missing, length(interactions), length(rel_zetaVals))
current_data   = Array{Union{Vector{Float64}, Missing}}(missing, length(interactions), length(rel_zetaVals))
vortex_rad_data = Array{Union{Float64, Missing}}(missing, length(interactions), length(rel_zetaVals))

# Storage arrays for fig 1 (density profiles)
r_xi_fig1        = Array{Union{Vector{Float64}, Missing}}(missing, length(interactions), length(rel_zetaVals))
rho_tot_norm_data = Array{Union{Vector{Float64}, Missing}}(missing, length(interactions), length(rel_zetaVals))
rho_s_norm_data  = Array{Union{Vector{Float64}, Missing}}(missing, length(interactions), length(rel_zetaVals))
rho_n_norm_data  = Array{Union{Vector{Float64}, Missing}}(missing, length(interactions), length(rel_zetaVals))
rho_imb_norm_data = Array{Union{Vector{Float64}, Missing}}(missing, length(interactions), length(rel_zetaVals))

if !read_from_table
    for j in eachindex(interactions), k in eachindex(rel_zetaVals)
        int = interactions[j]
        T = temperatures[1]
        β = (T == 0) ? 1000 : 1/(T*find_Tc(int))
        ζ = rel_zetaVals[k] * zetaCrit[j]
        println("int = $int, T/T_F = $(1/β), ζ = $ζ")

        (; Δ, μ, A_itp, Dt_itp, Ct, Et, Q, Rt, G, ξ2) = calculate_mf_eft_parameters(β,int,ζ)
        ξ2 = sqrt(2*Ct /(G*Δ^2))

        relaxed_data = readdlm(normpath(joinpath(PROJECT_ROOT,"data/relaxed_profile_$(bc)_int=$(int)_T=$(T)Tc_ζ=$(rel_zetaVals[k])ζc.csv")),',')
        u = Δ .* relaxed_data[:,2]
        r = relaxed_data[:,1]

        Δ0, μ0 = find_Δ_μ(1e3, int, 0.0)

        # Fig 6 quantities
        f       = u ./ Δ0
        current = (u.^2 / Δ0^2) ./ (r / ξ2)
        current[1] = 0.0  # u ~ r near the core, so j ~ r -> 0
        idx     = argmax(current[begin+1:end])

        r_xi_fig6[j,k]      = r/ξ2
        f_data[j,k]         = f
        current_data[j,k]   = current
        vortex_rad_data[j,k] = r[idx+1] / ξ2

        open("$PROJECT_ROOT/JLTP2026/tables/fig6_profiles_T=0_int=$(int)_zeta=$(rel_zetaVals[k])zetac.csv", "w") do io
            write(io, "# This data corresponds to (a subplot of) figure 6 of https://arxiv.org/abs/2512.22099\n# If using this data, please cite the above paper.\n")
            writedlm(io, permutedims(["r/xi", "f(r)=u/Delta0", "j(r)=(u/Delta0)^2/(r/xi)"]), ',')
            writedlm(io, hcat(r/ξ2, f, current), ',')
        end

        # Fig 1 quantities
        ρ_∞  = 1/6π^2
        ρ_tot = 1/2*broadcast(Φ->numberDensity_mf(β, μ, ζ, Φ)[1], u)
        ρ_s   = @. 2*Ct*u^2
        ρ_n   = ρ_tot - ρ_s
        Δρ    = 1/2*broadcast(Φ->imbalanceDensity_mf(β, μ, ζ, Φ)[1], u)

        r_xi_fig1[j,k]        = r/ξ2
        rho_tot_norm_data[j,k] = ρ_tot/ρ_∞
        rho_s_norm_data[j,k]  = ρ_s/ρ_∞
        rho_n_norm_data[j,k]  = ρ_n/ρ_∞
        rho_imb_norm_data[j,k] = Δρ/ρ_∞

        open("$PROJECT_ROOT/JLTP2026/tables/fig1_densities_T=0_int=$(int)_zeta=$(rel_zetaVals[k])zetac.csv", "w") do io
            write(io, "# This data corresponds to (a subplot of) figure 1 of https://arxiv.org/abs/2512.22099\n# If using this data, please cite the above paper.\n")
            writedlm(io, permutedims(["r/xi", "rho_tot/rho_inf", "rho_s/rho_inf", "rho_n/rho_inf", "rho_imbalance/rho_inf"]), ',')
            writedlm(io, hcat(r/ξ2, ρ_tot/ρ_∞, ρ_s/ρ_∞, ρ_n/ρ_∞, Δρ/ρ_∞), ',')
        end
    end
else
    for j in eachindex(interactions), k in eachindex(rel_zetaVals)
        int = interactions[j]

        tbl6 = readdlm("$PROJECT_ROOT/JLTP2026/tables/fig6_profiles_T=0_int=$(int)_zeta=$(rel_zetaVals[k])zetac.csv", ',', comments=true, skipstart=3)
        r_xi_fig6[j,k]      = Float64.(tbl6[:,1])
        f_data[j,k]         = Float64.(tbl6[:,2])
        current_data[j,k]   = Float64.(tbl6[:,3])
        idx = argmax(current_data[j,k][begin+1:end])
        vortex_rad_data[j,k] = r_xi_fig6[j,k][idx+1]

        tbl1 = readdlm("$PROJECT_ROOT/JLTP2026/tables/fig1_densities_T=0_int=$(int)_zeta=$(rel_zetaVals[k])zetac.csv", ',', comments=true, skipstart=3)
        r_xi_fig1[j,k]        = Float64.(tbl1[:,1])
        rho_tot_norm_data[j,k] = Float64.(tbl1[:,2])
        rho_s_norm_data[j,k]  = Float64.(tbl1[:,3])
        rho_n_norm_data[j,k]  = Float64.(tbl1[:,4])
        rho_imb_norm_data[j,k] = Float64.(tbl1[:,5])
    end
end

include("initialize_makie.jl")

# Fig 6: order parameter, current and vortex radius (T=0)

fig = Figure(size = (490pt,200pt), figure_padding = (0,2pt,1pt,0))

ylabels = [L"f(r)", L"\frac{j(r)}{\rho_{s,\infty} \hbar/(2m\xi)} "]

ax = [Axis(fig[i,j], xlabel=L"r/\xi", ylabel=ylabels[i]) for i in 1:2, j in eachindex(interactions)]

linewidths = [0.75pt, 1.125pt, 1.5pt, 0.5pt]

zeta_colors = [ColorSchemes.tab10[1],ColorSchemes.tab10[2],ColorSchemes.tab10[5],ColorSchemes.tab10[3]]

for j in eachindex(interactions)

    if j != 1
        hideydecorations!(ax[1,j], grid=false)
        linkyaxes!(ax[1,1], ax[1,j])
        hideydecorations!(ax[2,j], grid=false)
        linkyaxes!(ax[2,1], ax[2,j])
    end

    xlims!(ax[1,j], 0, 3)
    ylims!(ax[1,j], -0.05, 1.05)
    xlims!(ax[2,j], 0, 3)
    ylims!(ax[2,j], -0.05, 0.6)

    linkxaxes!(ax[1,j], ax[2,j])
    int = interactions[j]

    for k in eachindex(rel_zetaVals)
        lines!(ax[1,j], r_xi_fig6[j,k], f_data[j,k], label=L"\zeta/\zeta_c = %$(rel_zetaVals[k])", linewidth=linewidths[1], color = zeta_colors[k])
        hidexdecorations!(ax[1,j], grid=false)

        lines!(ax[2,j], r_xi_fig6[j,k], current_data[j,k], label=L"\zeta/\zeta_c = %$(rel_zetaVals[k])", linewidth=linewidths[1], color = zeta_colors[k])

        vlines!(ax[2,j], [vortex_rad_data[j,k]], linestyle=:dot, color = zeta_colors[k], alpha = 0.5)
    end

    # asymptote (ζ = 0)
    ζ = rel_zetaVals[1] * zetaCrit[j]
    β = 1e3
    (; Δ, μ, A_itp, Dt_itp, Ct, Et, Q, Rt, G, ξ2) = calculate_mf_eft_parameters(β,int,ζ)
    ξ2 = sqrt(2*Ct /(G*Δ^2))
    r_xi_asym = r_xi_fig6[j,1]

    lines!(ax[1,j], r_xi_asym, (@. (1-1/4*(1/r_xi_asym)^2)), label="asymptote", linewidth=linewidths[1], linestyle=:dash, color = :gray)
end

ax_radius = Axis(fig[1:2, length(interactions)+1], xlabel=L"ζ/ζ_c", ylabel=L"R_v/\xi", width = 110pt)

int_colors = mass_colors

data_radius = readdlm(normpath(joinpath(PROJECT_ROOT,"JLTP2026/tables/fig6b_vortex_radius_T=0.csv")),',',skipstart=3)
relzeta_radius = data_radius[:,1]
vortex_radii = data_radius[:,2:end]

@assert size(vortex_radii,2) == length(interactions)

for j in eachindex(interactions)
    scatter!(ax_radius, relzeta_radius, vortex_radii[:,j], label=L"(k_F a_s)^{-1} = %$(interactions[j])", color = int_colors[j])
end

axislegend(ax_radius, position = :lt, rowgap = 0, colgap=0, height=60pt, width = 90pt, margins = (0.5pt,0.5pt,0.5pt,0.5pt))

fig[0, 1:length(interactions)] = Legend(fig, ax[1,1], orientation = :horizontal)

rowgap!(fig.layout, 10)

for j in eachindex(interactions)
    Label(fig[1,j,Top()], L"(k_F a_s)^{-1} = %$(interactions[j])", halign = :center, padding = (0,0,4,0))
end

Label(fig[1,1,TopLeft()], "(a)")
Label(fig[1,length(interactions)+1,TopLeft()], "(b)")

save("$PROJECT_ROOT/JLTP2026/figures/fig6_profiles_and_current_T=0.pdf", fig)


# Fig 1: density profiles (T=0)

fig = Figure(size = (490pt,280pt), figure_padding = (0,2pt,2pt,1pt))

rmax = 4.0
zeta_indices = [1,2,4] # indices of rel_zetaVals to plot

ax = [k in zeta_indices ? Axis(fig[findall(x->x==k, zeta_indices)[1],j], xlabel=L"r/\xi", ylabel=L"\rho/\rho_{\text{tot},\infty}") : nothing for k in eachindex(rel_zetaVals), j in eachindex(interactions)]

linewidths = [0.75pt, 1.125pt, 1.5pt, 0.5pt]

for j in eachindex(interactions), k in zeta_indices
    if k ≠ length(rel_zetaVals)
        hidexdecorations!(ax[k,j], grid=false)
    end
    if j ≠ 1
        hideydecorations!(ax[k,j], grid=false)
    end

    int = interactions[j]
    ζ = rel_zetaVals[k] * zetaCrit[j]
    β = 1e3
    Δ, μ = find_Δ_μ(β, int, ζ)

    if Δ ≠ 0.0 && !ismissing(r_xi_fig1[j,k])
        (; Δ, μ, A_itp, Dt_itp, Ct, Et, Q, Rt, G, ξ2) = calculate_mf_eft_parameters(β,int,ζ)

        ρ_∞ = 1/6π^2
        ρ_s∞ = 2*Ct*Δ^2
        ρ_n∞ = ρ_∞ - ρ_s∞
        δρ_n∞ = ρ_s∞ - Δ^2 * Dt_int(β, μ, ζ, Δ)[1]

        r_xi        = r_xi_fig1[j,k]
        rho_tot_norm = rho_tot_norm_data[j,k]
        rho_s_norm  = rho_s_norm_data[j,k]
        rho_n_norm  = rho_n_norm_data[j,k]
        rho_imb_norm = rho_imb_norm_data[j,k]

        xlims!(ax[k,j], 0, rmax)
        ylims!(ax[k,j], -0.05, 1.05)

        lines!(ax[k,j], r_xi, rho_tot_norm, label="total", color=mass_colors[1], linewidth=linewidths[1])
        lines!(ax[k,j], r_xi, rho_s_norm, label="superfluid", color=mass_colors[2], linewidth=linewidths[2])
        lines!(ax[k,j], r_xi, rho_n_norm, label="normal", color=mass_colors[3], linewidth=linewidths[3])
        lines!(ax[k,j], r_xi, rho_imb_norm, label="imbalance", color = :black, linestyle = :dot)

        nanif(x,bool) = bool ? NaN : x

        lines!(ax[k,j], r_xi, nanif.(ρ_s∞/ ρ_∞ .* (1 .- 1/2 * r_xi.^-2), r_xi .<= 1), color = asymptotic_colors[2], linestyle=:dash, alpha = 0.8)
        lines!(ax[k,j], r_xi, nanif.(δρ_n∞ / ρ_∞ * 1/2 * r_xi.^-2, r_xi .<= 1), color = asymptotic_colors[3], linestyle=:dash, alpha =0.8)
    end
end

elem_t = LineElement(color = mass_colors[1], linewidth = linewidths[1])
elem_s = [LineElement(color = mass_colors[2], linewidth = linewidths[2], points = Point2f[(0, 0.6), (1, 0.6)]), LineElement(color = asymptotic_colors[2], points = Point2f[(0, 0.4), (1, 0.4)], linestyle=:dash, alpha = 0.8)]
elem_n = [LineElement(color = mass_colors[3], linewidth = linewidths[3], points = Point2f[(0, 0.6), (1, 0.6)]), LineElement(color = asymptotic_colors[3], points = Point2f[(0, 0.4), (1, 0.4)], linestyle=:dash, alpha = 0.8)]
elem_i = LineElement(color = :black, linewidth = linewidths[3], linestyle = :dot)
Legend(fig[0, :], [elem_t, elem_s, elem_n, elem_i], ["total", "superfluid", "normal", "imbalance"], orientation = :horizontal)

rowgap!(fig.layout, 10)

for k in zeta_indices, j in eachindex(interactions)
    text!(
        ax[k,j], 0.975, 0.14,
        text = L"ζ/ζ_c = %$(rel_zetaVals[k])",
        align = (:right, :bottom),
        space = :relative,
        fontsize = 9pt,
    )
end

for j in eachindex(interactions)
    Label(fig[1,j,Top()], L"(k_F a_s)^{-1} = %$(interactions[j])", halign = :center, padding = (0,0,4,0))
end

save("$PROJECT_ROOT/JLTP2026/figures/fig1_densities_T=0.pdf", fig)
