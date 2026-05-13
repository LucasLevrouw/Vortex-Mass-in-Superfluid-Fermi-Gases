### Code written by Lucas Levrouw (2026)
### This script generates figure 8 of https://arxiv.org/abs/2512.22099
### If using this code, please cite the above paper.

using VortexMass; const PROJECT_ROOT = pkgdir(VortexMass)
using DelimitedFiles

# Create necessary directories
mkpath("$PROJECT_ROOT/JLTP2026/tables/")
mkpath("$PROJECT_ROOT/JLTP2026/figures/")

read_from_table = false  # set to true to read from saved tables instead of computing

zeta_crit_data = readdlm(normpath(joinpath(PROJECT_ROOT, "JLTP2026/critical_ζ.dat")), '\t')
zeta_crit_interactions = zeta_crit_data[:,1]
zeta_crit_vals = zeta_crit_data[:,2]


temperatures = 0.0:0.1:0.9
interactions = [-1.0, 0.0, 1.0]
rel_zetaVals = 0.0:0.01:1.0

zetaCrit = Array{Union{Float64, Missing}}(missing, length(interactions))
critical_temperatures = find_Tc.(interactions)

for j in eachindex(interactions)
    int = interactions[j]
    index_zetacrit = findfirst(int .== zeta_crit_interactions)
    if index_zetacrit === nothing
        @warn "No ζ_c data for int = $int, code will not work."
        continue
    end
    zetaCrit[j] = zeta_crit_vals[index_zetacrit]
end

ρ_∞ = 1/6π^2
superfluid_densities = Array{Union{Float64, Missing}}(missing, length(temperatures), length(interactions), length(rel_zetaVals))
healing_lengths = copy(superfluid_densities)

if !read_from_table
    for k in eachindex(rel_zetaVals)
        println("Calculating for ζ/ζ_c = $(rel_zetaVals[k])")

        for j in eachindex(interactions), i in eachindex(temperatures)
            int = interactions[j]
            T = temperatures[i]
            Tc = critical_temperatures[j]
            β = (T == 0) ? 1000 : 1/(T*Tc)

            if rel_zetaVals[k] == 1.0
                superfluid_densities[i,j,k] = T == 0 ? ρ_∞ : 0.0
                continue
            end

            ζ = zetaCrit[j] * rel_zetaVals[k]
            Δ, μ = find_Δ_μ(β, int, ζ)

            if Δ == 0.0
                @warn "Δ == 0. Cannot calculate EFT coefficients in normal state."
                superfluid_densities[i,j,k] = 0.0
                continue
            end

            Ct = C_int(β,μ,ζ,Δ)[1]
            Dt = Dt_int(β, μ, ζ, Δ)[1]
            G  = G_int(β,μ,ζ,Δ)[1]

            ξ2 = sqrt(2Ct /(G*Δ^2))
            healing_lengths[i,j,k] = ξ2

            ρ_s_∞ = 2*Ct*Δ^2
            superfluid_densities[i,j,k] = ρ_s_∞
        end
    end

    for j in eachindex(interactions)
        int = interactions[j]
        open("$PROJECT_ROOT/JLTP2026/tables/fig8_sf_density_and_healing_length_int=$(int).csv", "w") do io
            write(io, "# This data corresponds to figure 8 of https://arxiv.org/abs/2512.22099\n# If using this data, please cite the above paper.\n")
            header = vcat(["zeta/zetac"],
                          ["rho_s/rho_inf_T=$(T)Tc" for T in temperatures],
                          ["xi_T=$(T)Tc"            for T in temperatures])
            writedlm(io, permutedims(header), ',')
            data = hcat(collect(rel_zetaVals),
                        [superfluid_densities[i,j,:] for i in eachindex(temperatures)]...,
                        [healing_lengths[i,j,:]      for i in eachindex(temperatures)]...)
            writedlm(io, data, ',')
        end
    end
else
    for j in eachindex(interactions)
        int = interactions[j]
        tbl = readdlm("$PROJECT_ROOT/JLTP2026/tables/fig8_sf_density_and_healing_length_int=$(int).csv", ',', Any, comments=true, skipstart=3)
        parse_cell(x) = x == "missing" ? missing : Float64(x)
        nT = length(temperatures)
        # columns: 1=zeta/zetac, 2..nT+1=rho_s, nT+2..2nT+1=xi
        for i in eachindex(temperatures)
            superfluid_densities[i,j,:] = parse_cell.(tbl[:, i+1])
            healing_lengths[i,j,:]      = parse_cell.(tbl[:, nT+i+1])
        end
    end
end

println("Plotting...")

include("initialize_makie.jl")
using MakieExtra

colors1 = ColorSchemes.tab10
colors2 = Makie.wong_colors()


# Plot all asymptotic values

temperature_indices = [1,3,5,7,9,10]

fig = Figure(size = (490pt,230pt), figure_padding = (0.5pt,5pt,2pt,1.5pt))

axes = Matrix{Axis}(undef, 3, length(interactions))
for j in eachindex(interactions)
    axes[1,j] = Axis(fig[1, j], xlabel = L"\zeta/\zeta_c", ylabel = L"\rho_{s,\infty}/\rho_{\text{tot},\infty}", limits = ((0.0,1.0),(nothing,nothing)), xticks = 0.0:0.2:1.0)
    axes[2,j] = Axis(fig[2, j], xlabel = L"\zeta/\zeta_c", ylabel = L"k_F \xi", limits = ((0.0,1.0),(1,30)), xticks = 0.0:0.2:1.0, yscale=log10,
    yticks = BaseMulTicks([1,3]), yminorticks=BaseMulTicks(2:9), yminorticksvisible=true, yminorticksize=1.5)

    hidexdecorations!(axes[1,j], grid=false)
    if j != 1
        hideydecorations!.(axes[1,j], grid=false)
        hideydecorations!.(axes[2,j], grid=false)
    end
end

linkaxes!(axes[1,:]...)
linkaxes!(axes[2,:]...)

for j in eachindex(interactions), i in temperature_indices
    lines!(axes[1,j], rel_zetaVals, superfluid_densities[i,j,:]./ρ_∞, label=L"T/T_c = %$(temperatures[i])", color = colors1[i])
    lines!(axes[2,j], rel_zetaVals, healing_lengths[i,j,:], label=L"T/T_c = %$(temperatures[i])", color = colors1[i])
end

for j in eachindex(interactions)
    Label(fig[1,j,Top()], L"(k_F a_s)^{-1} = %$(Int64(interactions[j]))", halign = :center, padding = (0,0,4,0))
end

Legend(fig[0,1:length(interactions)], axes[1], orientation = :horizontal)

save("$(PROJECT_ROOT)/JLTP2026/figures/fig8_sf_densities_and_healing_lengths_as_a_function_of_zeta.pdf", fig)
