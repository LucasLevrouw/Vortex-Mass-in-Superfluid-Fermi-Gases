### Code written by Lucas Levrouw (2026)
### This script generates figure 6b of https://arxiv.org/abs/2512.22099
### If using this code, please cite the above paper.


# Run calculate-profiles.jl first to generate the data

using Revise
using VortexMass; const PROJECT_ROOT = pkgdir(VortexMass)
using DelimitedFiles

# Create necessary directories
mkpath("$PROJECT_ROOT/JLTP2026/tables/")
mkpath("$PROJECT_ROOT/JLTP2026/figures/")

read_from_table = false  # set to true to read vortex radii from saved table

zeta_crit_data = readdlm(normpath(joinpath(PROJECT_ROOT, "JLTP2026/critical_ζ.dat")), '\t')
zeta_crit_interactions = zeta_crit_data[:,1]
zeta_crit_vals = zeta_crit_data[:,2]


interactions = [-1.0, 0.0, 1.0]
rel_zetaVals = 0.0:0.05:0.95
zetaCrit = zeros(length(interactions))

bc = "neumann"

core_radii = Array{Union{Float64, Missing}}(missing, length(interactions), length(rel_zetaVals))
healing_lengths = copy(core_radii)
superfluid_densities = copy(core_radii)

if !read_from_table
    # Write heading of table
    open("$PROJECT_ROOT/JLTP2026/tables/fig6b_vortex_radius_T=0.csv", "w") do io
        write(io, "# This data corresponds to (an inset of) figure 6b of https://arxiv.org/abs/2512.22099\n# If using this data, please cite the above paper.\n")
        writedlm(io, hcat(["zeta/zeta_c"], permutedims(["R_v/xi [(k_F a_s)^-1=$(int)]" for int in interactions])), ',')
    end

    for k in eachindex(rel_zetaVals)

        for j in eachindex(interactions)
        int = interactions[j]
        T = 0.0
        β = 1000.0


        index_zetacrit = findfirst(int .== zeta_crit_interactions)

        if index_zetacrit === nothing
            println("No ζ_c data for int = $int, skipping...")
            continue
        end
        zeta_crit = zeta_crit_vals[index_zetacrit]

        zetaCrit[j] = zeta_crit

        ζ = rel_zetaVals[k] * zeta_crit

        (; Δ, μ, A_itp, Dt_itp, Ct, Et, Q, Rt, G, ξ2) = calculate_mf_eft_parameters(β,int,ζ)

        data_path = normpath(joinpath(PROJECT_ROOT,"data/relaxed_profile_$(bc)_int=$(int)_T=$(T)Tc_ζ=$(rel_zetaVals[k])ζc.csv"))
        if !isfile(data_path)
            println("Data file not found for int = $(int), T/Tc = $(T), ζ/ζ_c = $(rel_zetaVals[k]), skipping...")
            continue
        end
        relaxed_data = readdlm(data_path,',')

        r = relaxed_data[:,1]
        u = Δ*relaxed_data[:,2]

        idx_max_current = argmax((2*Ct*u.^2 ./ r)[begin+1:end])

        core_radii[j,k] = r[idx_max_current]
        healing_lengths[j,k] = ξ2

        end
        # Write data to table
        open("$PROJECT_ROOT/JLTP2026/tables/fig6b_vortex_radius_T=0.csv", "a") do io
            writedlm(io, hcat(rel_zetaVals[k], permutedims(core_radii[:,k] ./ healing_lengths[:,k])), ',')
        end
    end
else
    data = readdlm("$PROJECT_ROOT/JLTP2026/tables/fig6b_vortex_radius_T=0.csv", ',', comments=true, skipstart=3)
    for k in eachindex(rel_zetaVals), j in eachindex(interactions)
        core_radii[j, k] = data[k, j+1]
        healing_lengths[j, k] = 1.0
    end
end

 
## Create plots and tables

include("initialize_makie.jl")
light_mass_colors = ColorSchemes.tab20[2:2:6]



### Vortex radius as a function of ζ

fig = Figure(size = (200pt,150pt), figure_padding = (0.5pt,5pt,2pt,1.5pt))
ax = Axis(fig[1, 1], xlabel = L"\zeta/\zeta_c", ylabel = L"R_v/\xi", limits = ((-0.05,1.0),(nothing,nothing)), xticks = 0.0:0.2:1.0)

markers1 = [:circle, :rect, :utriangle]
markers2 = [:star5, :diamond, :dtriangle]

for j in eachindex(interactions)
    scatter!(ax, rel_zetaVals, core_radii[j,:]./healing_lengths[j,:], label = L"(k_F a_s)^{-1} = %$(interactions[j])")
end

axislegend(ax, position = :lt)

save("$(PROJECT_ROOT)/JLTP2026/figures/fig6b_vortex_radius_T=0.pdf", fig)

println("Done.")
