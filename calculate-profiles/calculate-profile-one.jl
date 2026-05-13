println("Starting script.")
println("Arguments received: $ARGS")

using VortexMass; const PROJECT_ROOT = pkgdir(VortexMass)
using OrdinaryDiffEq
using Optim
using DelimitedFiles

# include(joinpath(PROJECT_ROOT, "JLTP2026/eft_mf_imbalance.jl"))


zeta_crit_data = readdlm(normpath(joinpath(PROJECT_ROOT, "JLTP2026/critical_ζ.dat")), '\t')
zeta_crit_interactions = zeta_crit_data[:,1]
zeta_crit_vals = zeta_crit_data[:,2]


R_ξ2 = 150
N = 7500
l = 1
bc = :neumann


T = parse(Float64,ARGS[1])
int = parse(Float64, ARGS[2])
rel_zeta = parse(Float64, ARGS[3])
if !(0 ≤ T ≤ 1 && 0 ≤ rel_zeta ≤ 1)
    error("Invalid arguments: T and rel_zeta must be in [0,1].")
end

datapath =  "$PROJECT_ROOT/data/"
filepath         = "$datapath/relaxed_profile_$(bc)_int=$(int)_T=$(T)Tc_ζ=$(rel_zeta)ζc.csv"
filepath_fenergy = "$datapath/free_energy_$(bc)_int=$(int)_T=$(T)Tc_ζ=$(rel_zeta)ζc.csv"

if isfile(filepath) && isfile(filepath_fenergy)
    println("Profile and free energy already calculated, skipping...")
    exit()
end

index_zetacrit = findfirst(int .== zeta_crit_interactions)

if index_zetacrit === nothing
    println("No ζ_c data for int = $int, skipping...")
    exit()
end

ζ = rel_zeta * zeta_crit_vals[index_zetacrit]
Tc = find_Tc(int) #critical_temperatures[j]

β = (T == 0) ? 1000 : 1/(T*Tc)

println("int = $int, T/T_F = $(1/β), ζ = $ζ")


mf_eft_params = calculate_mf_eft_parameters(β,int,ζ)
ξ2 = mf_eft_params.ξ2

if mf_eft_params.Δ == 0.0
    println("No superfluid solution found for int = $int, T/T_F = $(1/β), ζ = $ζ. Stopping script.")
    exit()
end

R = R_ξ2 * ξ2

t = 2e4

println("Evolving in imaginary time for tau = $t\nStarting computation...")

callback, saved_values = free_energy_callback(mf_eft_params)
# If not converged, try AutoTsit5(TRBDF2())
@time profile = calculate_profile(mf_eft_params, R, N, t, bc; callback, solver = AutoTsit5(Rodas5())) 

mkpath(datapath)
writedlm(filepath, profile, ',')
writedlm(filepath_fenergy, [saved_values.t saved_values.saveval], ',')
