### Code written by Lucas Levrouw (2025)
### This script generates figures 2 and 4 of https://arxiv.org/abs/2505.12590
### If using this code, please cite the above paper.



# Run calculate-profiles.jl first to generate the data

using Revise
using VortexMass; const PROJECT_ROOT = pkgdir(VortexMass)
using OrdinaryDiffEq
using DelimitedFiles
using Trapz
using LaTeXStrings

# Create necessary directories
mkpath("$PROJECT_ROOT/PRA2025/tables/")
mkpath("$PROJECT_ROOT/PRA2025/figures/")

read_from_table = false  # set to true to read masses from saved tables

temperatures = [0.0, 0.5, 0.9]
interactions = -2.0:0.2:3.0
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

 
critical_temperatures = find_Tc.(interactions)

L_kF = 150
L0_kF = 15


associated_masses = zeros(length(temperatures), length(interactions))
internal_masses = copy(associated_masses)
associated_masses0 = copy(associated_masses)
internal_masses0 = copy(associated_masses)

mass_scales = copy(associated_masses)
healing_lengths = copy(associated_masses)
superfluid_densities = copy(associated_masses)
coefficients_internal = copy(associated_masses)

correction_factor_a = copy(associated_masses)
correction_factor_i = copy(associated_masses)

sound_speed = copy(associated_masses)

using NumericalIntegration



for j in eachindex(interactions), i in eachindex(temperatures)
    int = interactions[j]
    T = temperatures[i]
    Tc = find_Tc(int)
    β = (T == 0) ? 1000 : 1/(T*Tc)
    ζ = 0.0
    (; Δ, μ, A_itp, Dt_itp, Ct, Et, Q, Rt, G, ξ2) = calculate_mf_eft_parameters(β,int,ζ)
    u = Δ*relaxed_profiles[i,j]

    healing_lengths[i,j] = ξ2
    local L = L_kF/ξ2
    local L0 = L0_kF/ξ2

    
    sound_speed[i,j] = sqrt(2Ct*G *Δ^2 /(Dt_itp(Δ^2)^2+2G*Q *Δ^2))

    r = r_array[i,j]

    ρ_tot = 1/2*broadcast(Φ->numberDensity_mf(β, μ, ζ, Φ)[1],u)
    ρ_s = @. 2*Ct*u^2
    ρ_n = ρ_tot - ρ_s
    
    ρ_∞ = 1/6π^2
    ρ_s_∞ = 2*Ct*Δ^2
    ρ_n_∞ = ρ_∞ - ρ_s_∞

    coeff_a = ρ_s_∞
    coeff_i = -Δ^2 * Dt_int(β, μ, ζ, Δ)[1] + ρ_s_∞ 

    superfluid_densities[i,j] = ρ_s_∞
    coefficients_internal[i,j] = coeff_i

    idx = findfirst(x->x≥L*ξ2, r)

    m_a = 2π * integrate(r[1:idx], (r.*(ρ_s_∞  .-ρ_s))[1:idx])
    m_i = 2π * integrate(r[1:idx], (r.*(ρ_n .- ρ_n_∞))[1:idx])

    associated_masses[i,j] = m_a
    internal_masses[i,j] = m_i

    
    idx0 = findfirst(x->x≥L0*ξ2, r)

    associated_masses0[i,j] =  2π * integrate(r[1:idx0], (r.*(ρ_s_∞  .-ρ_s))[1:idx0])# m_a + 2π*ξ1^2 * coeff_a*log(L0* ξ1/r[end]) 
    internal_masses0[i,j] =  2π * integrate(r[1:idx0], (r.*(ρ_n .- ρ_n_∞))[1:idx0]) # m_i + 2π*ξ1^2 * coeff_i*log(L0*ξ1/r[end])

    mass_scales[i,j] = π*ξ2^2 * ρ_s_∞

    m_a0 = 2π * integrate(r, (r.*(ρ_s_∞  .-ρ_s))) + π*ξ2^2 * coeff_a*log(ξ2/r[end])
    m_i0 = 2π * integrate(r, (r.*(ρ_n .- ρ_n_∞))) + π*ξ2^2 * coeff_i*log(ξ2/r[end])
 

    correction_factor_a[i,j] = exp(-m_a0/(π*ξ2^2 * coeff_a))
    correction_factor_i[i,j] = exp(-m_i0/(π*ξ2^2 * coeff_i))
end

total_masses = associated_masses .+ internal_masses
total_masses0 = associated_masses0 .+ internal_masses0

systemsizes = [10, 150, Inf]
R_labels = [L"k_F R=10", L"k_F R=150", L"k_F R=\infty"]
mass_ratios = zeros(length(interactions), length(systemsizes)) # only for T = 0

for (k, R) in enumerate(systemsizes)
    for j in eachindex(interactions)
        int = interactions[j]
        β = 1000
        ζ = 0.0
  
        (; Δ, μ, A_itp, Dt_itp, Ct, Et, Q, Rt, G, ξ2) = calculate_mf_eft_parameters(β,int,ζ)
        if Δ == 0.0
            println("Found  Δ = 0. Skipping int = $int, T = $T Tc, ζ = $ζc.")
            continue
        end
        u = Δ*relaxed_profiles[1,j]
        r = r_array[1,j]
    
        ρ_tot = 1/2*broadcast(Φ->numberDensity_mf(β, μ, ζ, Φ)[1],u)
        ρ_s = @. 2*Ct*u^2
        ρ_n = ρ_tot - ρ_s
        
        ρ_∞ = 1/6π^2
        ρ_s_∞ = 2*Ct*Δ^2
        ρ_n_∞ = ρ_∞ - ρ_s_∞

        coeff_a = 2*Ct*Δ^2
        coeff_i = -Δ^2 * Dt_int(β, μ, ζ, Δ)[1] + 2*Ct*Δ^2

        if R == Inf
            mass_ratios[j,k] = coeff_i/(coeff_a + coeff_i)
        else
            if R > r[end] && !(R ≈ r[end])
                @show R
                throw(ArgumentError("R must be smaller than the numerical domain size"))
            end
            idx = (R ≈ r[end]) ? length(r) : findfirst(x->x≥R, r)
            m_a = 2π * integrate(r[1:idx], (r.*(ρ_s_∞  .-ρ_s))[1:idx])
            m_i = 2π * integrate(r[1:idx], (r.*(ρ_n .- ρ_n_∞))[1:idx])

            mass_ratios[j,k] = m_i/(m_a + m_i)
        end

    end
end

mass_ratios_temp = zeros(length(temperatures), length(interactions))

for i in eachindex(temperatures)
    for j in eachindex(interactions)
        int = interactions[j]
        ζ = 0
        T = temperatures[i]
        Tc = find_Tc(int)
        β = (T == 0) ? 1000 : 1/(T*Tc)
  
        (; Δ, μ, A_itp, Dt_itp, Ct, Et, Q, Rt, G, ξ2) = calculate_mf_eft_parameters(β,int,ζ)
        if Δ == 0.0
            println("Found  Δ = 0. Skipping int = $int, T = $T Tc, ζ = $ζc.")
            continue
        end

        u = Δ*relaxed_profiles[1,j]
        r = r_array[1,j]
    
        ρ_tot = 1/2*broadcast(Φ->numberDensity_mf(β, μ, ζ, Φ)[1],u)
        ρ_s = @. 2*Ct*u^2
        ρ_n = ρ_tot - ρ_s
        
        ρ_∞ = 1/6π^2
        ρ_s_∞ = 2*Ct*Δ^2
        ρ_n_∞ = ρ_∞ - ρ_s_∞

        coeff_a = 2*Ct*Δ^2
        coeff_i = -Δ^2 * Dt_int(β, μ, ζ, Δ)[1] + 2*Ct*Δ^2

            mass_ratios_temp[i,j] = coeff_i/(coeff_a + coeff_i)

    end
end


systemsizes_cont = 1.0:1.0:20.0
interaction_indices = [6,11,16] # [1,6,11]
interactions_lim = interactions[interaction_indices]
correction_factor_a_sz = zeros(length(interaction_indices), length(systemsizes_cont))
correction_factor_i_sz = copy(correction_factor_a_sz)

for (k, R) in enumerate(systemsizes_cont)
    for (j, j_prev) in enumerate(interaction_indices)
        int = interactions[j_prev]
        β = 1000
        ζ = 0.0
  
        (; Δ, μ, A_itp, Dt_itp, Ct, Et, Q, Rt, G, ξ2) = calculate_mf_eft_parameters(β,int,ζ)
        u = Δ*relaxed_profiles[1,j_prev]
        r = r_array[1,j_prev]
    
        ρ_tot = 1/2*broadcast(Φ->numberDensity_mf(β, μ, ζ, Φ)[1],u)
        ρ_s = @. 2*Ct*u^2
        ρ_n = ρ_tot - ρ_s
        
        ρ_∞ = 1/6π^2
        ρ_s_∞ = 2*Ct*Δ^2
        ρ_n_∞ = ρ_∞ - ρ_s_∞

        coeff_a = 2*Ct*Δ^2
        coeff_i = -Δ^2 * Dt_int(β, μ, ζ, Δ)[1] + 2*Ct*Δ^2

        idx = findfirst(x->x≥R*ξ2, r)

        m_a0 = 2π * integrate(r[1:idx], (r.*(ρ_s_∞  .-ρ_s))[1:idx]) + π*ξ2^2 * coeff_a*log(ξ2/r[idx])
        m_i0 = 2π * integrate(r[1:idx], (r.*(ρ_n .- ρ_n_∞))[1:idx]) + π*ξ2^2 * coeff_i*log(ξ2/r[idx])
     
    
        correction_factor_a_sz[j,k] = exp(-m_a0/(π*ξ2^2 * coeff_a))
        correction_factor_i_sz[j,k] = exp(-m_i0/(π*ξ2^2 * coeff_i))
    end
end

else  # read_from_table
    ni = length(interactions)
    nT = length(temperatures)

    data2a = readdlm("$PROJECT_ROOT/PRA2025/tables/fig2a_vortex_mass_T=0.csv",       ',', comments=true, skipstart=3)
    data2b = readdlm("$PROJECT_ROOT/PRA2025/tables/fig2b_vortex_mass_ratio_T=0.csv", ',', comments=true, skipstart=3)
    data2c = readdlm("$PROJECT_ROOT/PRA2025/tables/fig2c_correction_factors_T=0.csv",',' , comments=true, skipstart=3)
    data4a = readdlm("$PROJECT_ROOT/PRA2025/tables/fig4a_vortex_mass_finite_temperature.csv",         ',', comments=true, skipstart=3)
    data4b = readdlm("$PROJECT_ROOT/PRA2025/tables/fig4b_vortex_mass_ratio_finite_temperature.csv",   ',', comments=true, skipstart=3)
    data4c = readdlm("$PROJECT_ROOT/PRA2025/tables/fig4c_correction_factors_finite_temperature.csv",  ',', comments=true, skipstart=3)

    associated_masses   = zeros(nT, ni); associated_masses[1, :] = data2a[:, 3] ./ 2
    internal_masses     = zeros(nT, ni); internal_masses[1, :]   = data2a[:, 4] ./ 2
    total_masses        = reshape(data4a[:, 2:end]', nT, ni) ./ 2
    total_masses0       = copy(total_masses)
    mass_ratios         = data2b[:, 2:end]
    mass_ratios_temp    = reshape(data4b[:, 2:end]', nT, ni)
    correction_factor_a = zeros(nT, ni)
    correction_factor_i = zeros(nT, ni)
    correction_factor_a[1, :] = data2c[:, 2]
    correction_factor_i[1, :] = data2c[:, 3]
    correction_factor_a[2, :] = data4c[:, 2]; correction_factor_a[3, :] = data4c[:, 4]
    correction_factor_i[2, :] = data4c[:, 3]; correction_factor_i[3, :] = data4c[:, 5]

    # Asymptotic quantities (fast: no profiles needed)
    healing_lengths        = zeros(nT, ni)
    superfluid_densities   = zeros(nT, ni)
    coefficients_internal  = zeros(nT, ni)
    for i in eachindex(temperatures), j in eachindex(interactions)
        int = interactions[j]; T = temperatures[i]
        β = (T == 0) ? 1000 : 1/(T*find_Tc(int))
        (; Δ, μ, Ct, ξ2) = calculate_mf_eft_parameters(β, int, 0.0)
        healing_lengths[i,j]       = ξ2
        superfluid_densities[i,j]  = 2*Ct*Δ^2
        coefficients_internal[i,j] = -Δ^2 * Dt_int(β, μ, 0.0, Δ)[1] + 2*Ct*Δ^2
    end

    systemsizes             = [10, 150, Inf]
    correction_factor_a_sz  = zeros(length(interaction_indices), length(systemsizes_cont))
    correction_factor_i_sz  = copy(correction_factor_a_sz)
end

## Create plots and tables

inch =  # in px
pt = 4/3 # in px (96 px == 1 inch)

using CairoMakie, LaTeXStrings

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
    Lines = (
        linewidth = 1.5,
    ),
    ScatterLines = (
        markersize = 7,
        linewidth = 1.5,
    ),
    )
set_theme!(merge(my_theme, theme_latexfonts()))
# default Plots.jl colors
using Colors, ColorSchemes
using Makie: wong_colors
plots_colors = [RGB{Float64}(0.0,0.6056031611752245,0.9786801175696073),
RGB{Float64}(0.8888735002725197,0.4356491903481899,0.2781229361419437),
RGB{Float64}(0.2422242978521988,0.6432750931576305,0.3044486515341153)]
colors = plots_colors
makie_colors = Makie.wong_colors()

asymptotic_colors = getindex.(Ref(ColorSchemes.tab10), 4:6)

temperature_colors = makie_colors
mass_colors = ColorSchemes.tab10
light_mass_colors = ColorSchemes.tab20[2:2:6]

R_colors = [ColorSchemes.tab10[7],ColorSchemes.tab10[9],ColorSchemes.tab10[8]]

ms = 5
lw = 1.5

fig = Figure(size = (510pt,200pt), figure_padding = (0,2,0,2))
ax1 = Axis(fig[1:2, 1], xlabel = L"(k_Fa_s)^{-1}", ylabel = L"\frac{M}{m \,k_F}", yscale = log10, width = 230pt)
ax1.xticks = -2:1:3

markers1 = [:circle, :rect, :utriangle]
markers2 = [:star5, :diamond, :dtriangle]

lines!(ax1, interactions, 2 * superfluid_densities[1,:] .* π .* healing_lengths[1,:].^2 , color = :gray)
lines!(ax1, interactions, 2 * coefficients_internal[1,:] .* π .* healing_lengths[1,:].^2 , color = :black, linestyle = :dash)

scatter!(ax1, interactions, 2 * total_masses[1, :], marker = :circle, color = mass_colors[1], label="total")

scatter!(ax1, interactions, 2 * associated_masses[1, :], marker = :rect, color = mass_colors[2], label="associated")

scatter!(ax1, interactions, 2 * internal_masses[1, :], marker = :utriangle, color = mass_colors[3], label="internal")

lines!(ax1, interactions, 2 * superfluid_densities[1,:] .* π .* healing_lengths[1,:].^2 .* log.(L_kF ./ (exp(-3/4) *  healing_lengths[1,:])), color = asymptotic_colors[2], linestyle = :dot)
lines!(ax1, interactions, 2 * coefficients_internal[1,:] .* π .* healing_lengths[1,:].^2 .* log.(L_kF ./ (exp(-3/4) *  healing_lengths[1,:])),  color = asymptotic_colors[3], linestyle = :dot)

# write corresponding table
if !read_from_table
    open("$PROJECT_ROOT/PRA2025/tables/fig2a_vortex_mass_T=0.csv", "w") do io
        write(io, "# This data corresponds to figure 2a of https://arxiv.org/abs/2505.12590\n# If using this data, please cite the above paper.\n")
        writedlm(io, permutedims(["(k_F a_s)^(-1)", "M_tot/(m k_F)", "M_a/(m k_F)", "M_i/(m k_F)", "M_a_ana/(m k_F)", "M_i_ana/(m k_F)", "pi xi^2 rho_s / (m k_F)", "pi xi^2 deltarho_n / (m k_F)"]), ',')
        writedlm(io, hcat(interactions, 2 * total_masses[1, :], 2 * associated_masses[1, :], 2 * internal_masses[1, :],
            2 * superfluid_densities[1,:] .* π .* healing_lengths[1,:].^2 .* log.(L_kF ./ (exp(-3/4) *  healing_lengths[1,:])),
            2 * coefficients_internal[1,:] .* π .* healing_lengths[1,:].^2 .* log.(L_kF ./ (exp(-3/4) *  healing_lengths[1,:])),
            2 * superfluid_densities[1,:] .* π .* healing_lengths[1,:].^2,
            2 * coefficients_internal[1,:] .* π .* healing_lengths[1,:].^2), ','
        )
    end
end

axislegend(ax1, rowgap = 0)

ax2 = Axis(fig[1, 2], xlabel = L"(k_Fa)^{-1}", ylabel = L"M_i/M_{tot}",)
ax2.xticks = -2:1:3
ax2.yticks = 0.0:0.25:0.5

markersR = [:utriangle, :dtriangle, nothing]

for k in eachindex(systemsizes)
    if systemsizes[k] == Inf
        lines!(ax2, interactions, mass_ratios[:,k], color = R_colors[k], label=R_labels[k])
    else
        scatter!(ax2, interactions, mass_ratios[:,k], color = R_colors[k], marker = markersR[k], label=R_labels[k])
    end
end
axislegend(ax2, position=:rt, rowgap = -5, height=60)

# write corresponding table
if !read_from_table
    open("$PROJECT_ROOT/PRA2025/tables/fig2b_vortex_mass_ratio_T=0.csv", "w") do io
        write(io, "# This data corresponds to figure 2b of https://arxiv.org/abs/2505.12590\n# If using this data, please cite the above paper.\n")
        writedlm(io, permutedims(["(k_F a_s)^(-1)", "M_i/M_tot (k_F R = 10)", "M_i/M_tot (k_F R = 150)", "M_i/M_tot (k_F R = infinity)"]), ',')
        writedlm(io, hcat(interactions, mass_ratios), ',')
    end
end

ax3 = Axis(fig[2, 2], xlabel = L"(k_Fa_s)^{-1}", ylabel = L"\alpha")
ax3.xticks = -2:1:3


scatter!(ax3, interactions, correction_factor_a[1, :], color = mass_colors[2], marker = :rect, label = L"\alpha_a")
scatter!(ax3, interactions, correction_factor_i[1, :], color = mass_colors[3], marker = :utriangle, label = L"\alpha_i")
lines!(ax3, [(-2,exp(-3/4)),(3,exp(-3/4))], color = :gray, linestyle = :dot) # !! Limits don't work

axislegend(ax3, position=:rc, rowgap = -5, height=50)

# write corresponding table
if !read_from_table
    open("$PROJECT_ROOT/PRA2025/tables/fig2c_correction_factors_T=0.csv", "w") do io
        write(io, "# This data corresponds to figure 2c of https://arxiv.org/abs/2505.12590\n# If using this data, please cite the above paper.\n")
        writedlm(io, permutedims(["(k_F a_s)^(-1)", "alpha_a", "alpha_i"]), ',')
        writedlm(io, hcat(interactions, correction_factor_a[1, :], correction_factor_i[1, :]), ',')
    end
end

rowgap!(fig.layout, 10)

Label(fig[1,1:2,TopLeft()], "(a)")
Label(fig[1,2,TopLeft()], "(b)", halign=:left)
Label(fig[2,2,TopLeft()], "(c)", halign=:left)

save("$(PROJECT_ROOT)/PRA2025/figures/fig2_vortex_mass_T=0.png", fig, px_per_unit = 10)
save("$(PROJECT_ROOT)/PRA2025/figures/fig2_vortex_mass_T=0.eps", fig, px_per_unit = 10)

println("Done")


fig = Figure(size = (490pt,210pt))
ax1 = Axis(fig[1:2, 1], xlabel = L"(k_Fa_s)^{-1}", ylabel = L"\frac{M_{tot}}{m\, k_F}", yscale = log10, width = 180pt)
ax1.xticks = -2:1:3

markers1 = [:circle, :rect, :utriangle]
markers2 = [:star5, :diamond, :dtriangle]

for i in eachindex(temperatures)
    scatter!(ax1, interactions, 2*total_masses[i, :], marker = markers1[i],
    color = temperature_colors[i], label=latexstring("T/T_c = $(temperatures[i])"))
end
axislegend(ax1, rowgap = 0)

# write corresponding table
if !read_from_table
    open("$PROJECT_ROOT/PRA2025/tables/fig4a_vortex_mass_finite_temperature.csv", "w") do io
        write(io, "# This data corresponds to figure 4a of https://arxiv.org/abs/2505.12590\n# If using this data, please cite the above paper.\n")
        writedlm(io, permutedims(["(k_F a_s)^(-1)", "M_tot/(m k_F) [T/T_c = 0]", "M_tot/(m k_F) [T/T_c = 0.5]", "M_tot/(m k_F) [T/T_c = 0.9]"]), ',')
        writedlm(io, hcat(interactions, 2 .* total_masses'), ','
        )
    end
end

ax2 = Axis(fig[1, 2], xlabel = L"(k_Fa)^{-1}", ylabel = L"M_i/M_{tot}",)
ax2.xticks = -2:1:3
ax2.yticks = 0.0:0.25:0.5

linestyles = [:solid, :dash, :dashdot]

for i in eachindex(temperatures)
    lines!(ax2, interactions, mass_ratios_temp[i,:], color = temperature_colors[i], label = latexstring("T/T_c = $(temperatures[i])"),linestyle = linestyles[i])
end
Legend(fig[1,3], ax2, tellheight=true)

if !read_from_table
    open("$PROJECT_ROOT/PRA2025/tables/fig4b_vortex_mass_ratio_finite_temperature.csv", "w") do io
        write(io, "# This data corresponds to figure 4b of https://arxiv.org/abs/2505.12590\n# If using this data, please cite the above paper.\n")
        writedlm(io, permutedims(["(k_F a_s)^(-1)", "M_i/M_tot (T/T_c = 0)", "M_i/M_tot (T/T_c = 0.5)", "M_i/M_tot (T/T_c = 0.9)"]), ',')
        writedlm(io, hcat(interactions, mass_ratios_temp'), ',')
    end
end



ax3 = Axis(fig[2, 2], xlabel = L"(k_Fa_s)^{-1}", ylabel = L"\alpha")
ax3.xticks = -2:1:3
linkxaxes!(ax2, ax3)

scatter!(ax3, interactions, correction_factor_a[2, :], color = mass_colors[2], marker=:rect, label = L"\alpha_a(T/T_c = 0.5)")
scatter!(ax3, interactions, correction_factor_i[2, :], color = mass_colors[3], marker = :utriangle, label = L"\alpha_i(T/T_c = 0.5)")

scatter!(ax3, interactions, correction_factor_a[3, :], color = asymptotic_colors[2], marker=:rect, label = L"\alpha_a(T/T_c = 0.9)")
scatter!(ax3, interactions, correction_factor_i[3, :], color = asymptotic_colors[3], marker=:utriangle, label = L"\alpha_i(T/T_c = 0.9)")

lines!(ax3, [(-2,exp(-3/4)),(3,exp(-3/4))], color = :gray, linestyle = :dot)
Legend(fig[2,3], ax3, tellheight = true)

# write corresponding table
if !read_from_table
    open("$PROJECT_ROOT/PRA2025/tables/fig4c_correction_factors_finite_temperature.csv", "w") do io
        write(io, "# This data corresponds to figure 4c of https://arxiv.org/abs/2505.12590\n# If using this data, please cite the above paper.\n")
        writedlm(io, permutedims(["(k_F a_s)^(-1)", "alpha_a(T/T_c = 0.5)", "alpha_i(T/T_c = 0.5)", "alpha_a(T/T_c = 0.9)", "alpha_i(T/T_c = 0.9)"]), ',')
        writedlm(io, hcat(interactions, correction_factor_a[2, :], correction_factor_i[2, :], correction_factor_a[3, :], correction_factor_i[3, :]), ',')
    end
end

rowgap!(fig.layout, 10)

Label(fig[1,1:2,TopLeft()], "(a)")
Label(fig[1,2,TopLeft()], "(b)", halign=:left)
Label(fig[2,2,TopLeft()], "(c)", halign=:left)

save("$(PROJECT_ROOT)/PRA2025/figures/fig4_vortex_mass-finite_temperature.png", fig, px_per_unit = 10)
save("$(PROJECT_ROOT)/PRA2025/figures/fig4_vortex_mass-finite_temperature.eps", fig, px_per_unit = 10)

println("Done")
