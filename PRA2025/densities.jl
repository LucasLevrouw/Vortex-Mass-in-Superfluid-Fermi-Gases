### Code written by Lucas Levrouw (2025)
### This script generates figures 1 and 3 of https://arxiv.org/abs/2505.12590
### If using this code, please cite the above paper.


# Run calculate-profiles.jl first to generate the data


using Revise
using VortexMass; const PROJECT_ROOT = pkgdir(VortexMass)
using OrdinaryDiffEq
using DelimitedFiles
using Trapz

# Create necessary directories
mkpath("$PROJECT_ROOT/PRA2025/tables/")
mkpath("$PROJECT_ROOT/PRA2025/figures/")

read_from_table = false  # set to true to read from saved tables instead of profiles

temperatures = [0.0, 0.5, 0.9]
interactions = -1.0:1.0:1.0

bc = "neumann"

if !read_from_table
    r_array  = Array{Union{Vector{Union{Float64, Missing}}, Missing}}(missing, length(temperatures), length(interactions))
    relaxed_profiles  = Array{Union{Vector{Union{Float64, Missing}}, Missing}}(missing, length(temperatures), length(interactions))

    for j in eachindex(interactions), i in eachindex(temperatures)
        int = interactions[j]
        T = temperatures[i]
        β = (T == 0) ? 1000 : 1/(T*find_Tc(int))
        (; Δ, μ, A_itp, Dt_itp, Ct, Et, Q, Rt, G, ξ2) = calculate_mf_eft_parameters(β,int,0.0)

        relaxed_data = readdlm(normpath(joinpath(PROJECT_ROOT,"data/relaxed_profile_$(bc)_int=$(int)_T=$(T)Tc_ζ=0.0ζc.csv")),',')
        r_array[i,j] = relaxed_data[:,1]
        relaxed_profiles[i,j] = relaxed_data[:,2]
    end
end
 
critical_temperatures = find_Tc.(interactions)

interactions_cont = -2.0:0.1:3.0

healing_length_2 =  zeros(length(temperatures),length(interactions_cont))
healing_length_bcs =  zeros(length(temperatures),length(interactions_cont))
healing_length_bec = fill(NaN, (length(temperatures),length(interactions_cont)))

for i in eachindex(temperatures), j in eachindex(interactions_cont)
    int = interactions_cont[j]
    T = temperatures[i]
    Tc = find_Tc(int)
    β = (T == 0) ? 1000 : 1/(T*Tc)
    ζ = 0.0
    (; Δ, μ, A_itp, Dt_itp, Ct, Et, Q, Rt, G, ξ2) = calculate_mf_eft_parameters(β,int,ζ)

    healing_length_2[i,j] = ξ2
    healing_length_bcs[i,j] = 2/(sqrt(3)) * 1/Δ 
    ρ_∞ = 1/6π^2
    g = 4π * 2*int^-1
    if g > 0
        healing_length_bec[i,j] = 1/sqrt(g*ρ_∞)
    end
end

L = 100 # box size in units of healing length

pt = 4/3 # in px (96 px == 1 inch)

using CairoMakie
using ColorSchemes
using Makie: wong_colors
my_theme = Theme(
    fontsize=8pt,
    labelsize=9pt,
    figure_padding = 0,
    Axis = (
        xticksize = 3,
        yticksize = 3,
    ),
        )
set_theme!(merge(my_theme, theme_latexfonts()))

temperature_colors = wong_colors()
mass_colors = ColorSchemes.tab10
asymptotic_colors = getindex.(Ref(ColorSchemes.tab10), 4:6) # ColorSchemes.tab20[2:2:end]
alpha = 1.0

rmax = 4


my_theme = Theme(
    fontsize=8pt,
    labelsize=9pt,
    figure_padding = 0,
    Axis = (
        xticksize = 3,
        yticksize = 3,
        xlabelsize = 9pt,
        ylabelsize = 9pt,
    ),
    Label = (
        fontsize = 9pt,
        ),
        )
set_theme!(merge(my_theme, theme_latexfonts()))

nanif(x,bool) = bool ? NaN : x


## FIGURE 1
## 1a. Total, superfluid and normal densities at T = 0
## 1b. Healing length at various temperatures

fig = Figure(size = (246pt, 420pt), figure_padding = (0,2,1,4)) # (left, right, bottom, top)

interactions = [-1, 0, 1]

ax = [Axis(fig[j,1:2], xlabel=L"r/\xi", ylabel=L"\rho/\rho_\infty") for j in eachindex(interactions)]
push!(ax, Axis(fig[4,1], xlabel=L"(k_F a_s)^{-1}", ylabel=L"k_F \xi", yscale=log10))


hidexdecorations!(ax[1], grid=false)
hidexdecorations!(ax[2], grid=false)


linewidths = [0.75pt, 1.125pt, 1.5pt]

for (j, int) in pairs(interactions)
    β = 1000
    ζ = 0.0
    (; Δ, μ, A_itp, Dt_itp, Ct, Et, Q, Rt, G, ξ2) = calculate_mf_eft_parameters(β,int,ζ)

    ρ_∞ = 1/6π^2
    ρ_s∞ = 2*Ct*Δ^2
    ρ_n∞ = ρ_∞ - ρ_s∞
    coeff_ind = 2*C_int(β, μ, ζ, Δ)[1]* Δ^2 # = ρ_s_∞
    coeff_int = -Δ^2 * Dt_int(β, μ, ζ, Δ)[1] + 2*C_int(β, μ, ζ, Δ)[1]* Δ^2
    nanif(x,bool) = bool ? NaN : x

    if !read_from_table
        u = Δ.*relaxed_profiles[1,j]
        r = r_array[1,j]
        r_xi = r/ξ2
        rho_tot_norm = 1/2*broadcast(Φ->numberDensity_mf(β, μ, ζ, Φ)[1],u) / ρ_∞
        rho_s_norm = (@. 2*Ct*u^2) / ρ_∞
        rho_n_norm = rho_tot_norm - rho_s_norm
    else
        tbl = readdlm("$PROJECT_ROOT/PRA2025/tables/fig1a_densities_T=0_(k_F a_s)^(-1)=$int.csv", ',', comments=true, skipstart=3)
        r_xi = Float64.(tbl[:,1])
        rho_tot_norm = Float64.(tbl[:,2])
        rho_s_norm = Float64.(tbl[:,3])
        rho_n_norm = Float64.(tbl[:,4])
    end

    xlims!(ax[j], 0, rmax)
    ylims!(ax[j], -0.05, 1.05)

    lines!(ax[j], r_xi, rho_tot_norm, label="total", color=mass_colors[1], linewidth=linewidths[1])
    lines!(ax[j], r_xi, rho_s_norm, label="superfluid", color=mass_colors[2], linewidth=linewidths[2])
    lines!(ax[j], r_xi, rho_n_norm, label="normal", color=mass_colors[3], linewidth=linewidths[3])

    lines!(ax[j], r_xi, nanif.(ρ_s∞/ ρ_∞ .* (1 .- 1/2 * r_xi.^-2), r_xi .< 1), color = asymptotic_colors[2], linestyle=:dash)
    lines!(ax[j], r_xi, nanif.(ρ_s∞/ ρ_∞ * (1/2) * r_xi.^2, r_xi .>= 1), color = asymptotic_colors[2], linestyle=:dashdot)

    lines!(ax[j], r_xi, nanif.(coeff_int / ρ_∞ * 1/2 * r_xi.^-2, r_xi .< 1), color = asymptotic_colors[3], linestyle=:dash)
    lines!(ax[j], r_xi, nanif.(coeff_int / ρ_∞ *(1 .- (1/2) * r_xi.^2), r_xi .>= 1), color = asymptotic_colors[3], linestyle=:dashdot)

    zeroif(x,bool) = bool ? 0.0 : x

    # Write corresponding tables
    if !read_from_table
        ρ_s_ana = zeroif.(ρ_s∞ .* (1 .- 1/2 * r_xi.^-2), r_xi .< 1) .+ zeroif.(ρ_s∞/ ρ_∞ * (1/2) * r_xi.^2, r_xi .>= 1)
        ρ_n_ana = zeroif.(coeff_int* 1/2 * r_xi.^-2, r_xi .< 1) + zeroif.(coeff_int / ρ_∞ *(1 .- (1/2) * r_xi.^2), r_xi .>= 1)
        open("$PROJECT_ROOT/PRA2025/tables/fig1a_densities_T=0_(k_F a_s)^(-1)=$int.csv", "w") do io
            write(io, "# This data corresponds to (a subplot of) figure 1a of https://arxiv.org/abs/2505.12590\n# If using this data, please cite the above paper.\n")
            writedlm(io, permutedims(["r/xi", "rho_tot/rho_inf", "rho_s/rho_inf", "rho_n/rho_inf", "rho_s_ana/rho_inf", "rho_n_ana/rho_inf"]), ',')
            writedlm(io, hcat(r_xi, rho_tot_norm, rho_s_norm, rho_n_norm, ρ_s_ana/ρ_∞, ρ_n_ana/ρ_∞ ), ',')
        end
    end

end

for j in eachindex(temperatures)
    lines!(ax[4],interactions_cont, healing_length_2[j,:], label=L"T/T_c = %$(temperatures[j])", color = temperature_colors[j])
end
lines!(ax[4],interactions_cont, healing_length_bcs[1,:], label="BCS", color = ColorSchemes.tab10[7], linestyle=:dash)
lines!(ax[4],interactions_cont,  healing_length_bec[1,:], label="BEC", color = ColorSchemes.tab10[9], linestyle=:dashdot)

    # Write corresponding table
    if !read_from_table
        open("$PROJECT_ROOT/PRA2025/tables/fig1b_healing_length.csv", "w") do io
            write(io, "# This data corresponds to figure 1b of https://arxiv.org/abs/2505.12590\n# If using this data, please cite the above paper.\n")
            writedlm(io, permutedims(["(k_F a_s)^(-1)", "k_F xi(T=0)", "k_F xi(T=0.5 T_c)", "k_F xi(T=0.9 T_c)", "k_F xi_BCS", "k_F xi_BEC"]), ',')
            writedlm(io, hcat(interactions_cont, healing_length_2[1,:], healing_length_2[2,:], healing_length_2[3,:], healing_length_bcs[1,:], healing_length_bec[1,:]), ',')
        end
    end

Makie.tightlimits!(ax[4])
Legend(fig[4, 2], ax[4], rowgap = 0, tellheight=true)


elem_t = LineElement(color = mass_colors[1], linewidth = linewidths[1])
elem_s = [LineElement(color = mass_colors[2], linewidth = linewidths[2], points = Point2f[(0, 0.6), (1, 0.6)]), LineElement(color = asymptotic_colors[2], points = Point2f[(0, 0.4), (1, 0.4)], linestyle=:dash)]
elem_n = [LineElement(color = mass_colors[3], linewidth = linewidths[3], points = Point2f[(0, 0.6), (1, 0.6)]), LineElement(color = asymptotic_colors[3], points = Point2f[(0, 0.4), (1, 0.4)], linestyle=:dash)]
Legend(fig[0, :], [elem_t, elem_s, elem_n], ["total", "superfluid", "normal"], orientation = :horizontal)


rowgap!(fig.layout, 10)

for (ax, label) in zip(ax[1:3], [L"(k_F a_s)^{-1} = -1", L"(k_F a_s)^{-1} = 0", L" (k_F a_s)^{-1} = 1"])
    text!(
        ax, 0.98, 0.12,
        text = label,
        align = (:right, :bottom),
        space = :relative,
        fontsize = 9pt,
    )
end

Label(fig[1,1,TopLeft()], "(a)", halign = :left)
Label(fig[4,1,TopLeft()], "(b)", halign = :left)

save("$PROJECT_ROOT/PRA2025/figures/fig1_densities_T=0_healing_length.png", fig, px_per_unit = 10)
save("$PROJECT_ROOT/PRA2025/figures/fig1_densities_T=0_healing_length.eps", fig, px_per_unit = 10)



## FIGURE 2: Total, superfluid and normal densities at various temperatures

fig = Figure(size = (246pt,246pt), figure_padding = (0,2,1,4)) # (left, right, bottom, top)



ax = [Axis(fig[j,i], xlabel=L"r/\xi", ylabel=L"\rho/\rho_\infty") for j in eachindex(interactions), i in 1:2]

hideydecorations!.(ax[:,2], grid=false)
hidexdecorations!.(ax[1:2,:], grid=false)

for i in 1:2, j in eachindex(interactions)
    int = interactions[j]

    T = temperatures[i+1]
    Tc = critical_temperatures[j]
    β = (T == 0) ? 1000 : 1/(T*Tc)

    ζ = 0.0
    (; Δ, μ, A_itp, Dt_itp, Ct, Et, Q, Rt, G, ξ2) = calculate_mf_eft_parameters(β,int,ζ)

    ρ_∞ = 1/6π^2
    ρ_s∞ = 2*Ct*Δ^2
    ρ_n∞ = ρ_∞ - ρ_s∞
    coeff_ind = 2*C_int(β, μ, ζ, Δ)[1]* Δ^2 # = ρ_s_∞
    coeff_int = -Δ^2 * Dt_int(β, μ, ζ, Δ)[1] + 2*C_int(β, μ, ζ, Δ)[1]* Δ^2

    if !read_from_table
        u = Δ.*relaxed_profiles[1,j]
        r = r_array[i,j]
        r_xi = r/ξ2
        rho_tot_norm = 1/2*broadcast(Φ->numberDensity_mf(β, μ, ζ, Φ)[1],u) / ρ_∞
        rho_s_norm = (@. 2*Ct*u^2) / ρ_∞
        rho_n_norm = rho_tot_norm - rho_s_norm
    else
        tbl = readdlm("$PROJECT_ROOT/PRA2025/tables/fig3_densities_T=$(T)_(k_F a_s)^(-1)=$int.csv", ',', comments=true, skipstart=3)
        r_xi = Float64.(tbl[:,1])
        rho_tot_norm = Float64.(tbl[:,2])
        rho_s_norm = Float64.(tbl[:,3])
        rho_n_norm = Float64.(tbl[:,4])
    end

    xlims!(ax[j, i], 0, rmax)
    ylims!(ax[j, i], -0.05, 1.05)

    lines!(ax[j, i], r_xi, rho_tot_norm, label="total", color=mass_colors[1], linewidth=linewidths[1])
    lines!(ax[j, i], r_xi, rho_s_norm, label="superfluid", color=mass_colors[2], linewidth=linewidths[2])
    lines!(ax[j, i], r_xi, rho_n_norm, label="normal", color=mass_colors[3], linewidth=linewidths[3])

    lines!(ax[j, i], r_xi, nanif.(ρ_s∞/ ρ_∞ .* (1 .-  1/2 * r_xi.^-2), r_xi .<= 1), color = asymptotic_colors[2], linestyle=:dash)
    lines!(ax[j, i], r_xi, nanif.(ρ_n∞/ρ_∞ .+ coeff_int/ρ_∞ * 1/2 * r_xi.^-2, r_xi .<= 1), color = asymptotic_colors[3], linestyle=:dash)

    lines!(ax[j, i], r_xi, nanif.(ρ_s∞/ ρ_∞ * (1/2) * r_xi.^2, r_xi .>= 1), color = asymptotic_colors[2], linestyle=:dashdot)
    lines!(ax[j, i], r_xi, ρ_n∞/ ρ_∞ .+ nanif.(coeff_int / ρ_∞ *(1 .- (1/2) * r_xi.^2), r_xi .>= 1), color = color = asymptotic_colors[3], linestyle=:dashdot)

    zeroif(x,bool) = bool ? 0.0 : x

    # Write corresponding tables
    if !read_from_table
        ρ_s_ana = zeroif.(ρ_s∞ .* (1 .- 1/2 * r_xi.^-2), r_xi .< 1) .+ zeroif.(ρ_s∞/ ρ_∞ * (1/2) * r_xi.^2, r_xi .>= 1)
        ρ_n_ana = zeroif.(ρ_n∞/ρ_∞ .+ coeff_int/ρ_∞ * 1/2 * r_xi.^-2, r_xi .<= 1) + zeroif.(ρ_n∞/ρ_∞ .+ coeff_int/ρ_∞ * 1/2 * r_xi.^-2, r_xi .<= 1)
        open("$PROJECT_ROOT/PRA2025/tables/fig3_densities_T=$(T)_(k_F a_s)^(-1)=$int.csv", "w") do io
            write(io, "# This data corresponds to (a subplot of) figure 3 of https://arxiv.org/abs/2505.12590\n# If using this data, please cite the above paper.\n")
            writedlm(io, permutedims(["r/xi", "rho_tot/rho_inf", "rho_s/rho_inf", "rho_n/rho_inf, rho_s_ana/rho_inf, rho_n_ana/rho_inf"]), ',')
            writedlm(io, hcat(r_xi, rho_tot_norm, rho_s_norm, rho_n_norm, ρ_s_ana/ρ_∞, ρ_n_ana/ρ_∞  ), ',')
        end
    end
end



elem_t = LineElement(color = mass_colors[1], linewidth = linewidths[1])
elem_s = [LineElement(color = mass_colors[2], linewidth = linewidths[2], points = Point2f[(0, 0.6), (1, 0.6)]), LineElement(color = asymptotic_colors[2], points = Point2f[(0, 0.4), (1, 0.4)], linestyle=:dash)]
elem_n = [LineElement(color = mass_colors[3], linewidth = linewidths[3], points = Point2f[(0, 0.6), (1, 0.6)]), LineElement(color = asymptotic_colors[3], points = Point2f[(0, 0.4), (1, 0.4)], linestyle=:dash)]
Legend(fig[0, :], [elem_t, elem_s, elem_n], ["total", "superfluid", "normal"], orientation = :horizontal)



rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

for (ax, label) in zip(ax[1:3, 1], [L"(k_F a_s)^{-1} = -1", L"(k_F a_s)^{-1} = 0", L" (k_F a_s)^{-1} = 1"])
    text!(
        ax, 0.98, 0.325,
        text = label,
        align = (:right, :bottom),
        space = :relative,
        fontsize = 9pt,
    )
    text!(
        ax, 0.98, 0.7,
        text = L"T/T_c = 0.5",
        align = (:right, :top),
        space = :relative,
        fontsize = 9pt,
    )
end

for (ax, label) in zip(ax[1:3, 2], [L"(k_F a_s)^{-1} = -1", L"(k_F a_s)^{-1} = 0", L" (k_F a_s)^{-1} = 1"])
    text!(
        ax, 0.98, 0.325,
        text = label,
        align = (:right, :bottom),
        space = :relative,
        fontsize = 9pt,
    )
    text!(
        ax, 0.98, 0.7,
        text = L"T/T_c = 0.9",
        align = (:right, :top),
        space = :relative,
        fontsize = 9pt,
    )
end



save("$PROJECT_ROOT/PRA2025/figures/fig3_densities_finite_temperature.png", fig, px_per_unit = 10)
save("$PROJECT_ROOT/PRA2025/figures/fig3_densities_finite_temperature.eps", fig, px_per_unit = 10)



