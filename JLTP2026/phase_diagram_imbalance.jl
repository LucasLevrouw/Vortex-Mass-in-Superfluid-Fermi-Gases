### Code written by Lucas Levrouw (2026)
### This script generates figure 5 of https://arxiv.org/abs/2512.22099
### It also generates critical_ζ.dat, which is used for all calculations with imbalance.
### If using this code, please cite the above paper.

using VortexMass; const PROJECT_ROOT = pkgdir(VortexMass)

using NLsolve, NonlinearSolve

betaVal = 1e4

intVals = -2.0:0.01:3.2 

 zetaVals =  [0.0,
 0.1,
 0.2,
 0.5,
 1.0,
 2.0,
 5.0,
10.0,
 ]

DeltaVals = zeros(length(intVals), length(zetaVals))
muVals = zeros(length(intVals), length(zetaVals))

for i = 1:length(intVals)
    println("int = $(intVals[i])")
    @time for k = 1:length(zetaVals)
        DeltaVals[i,k], muVals[i,k] = find_Δ_μ(betaVal,intVals[i],zetaVals[k])#, LevenbergMarquardt())
    end
end


 # Function to find critical interaction strength

 using NonlinearSolve, SciMLBase
 
 function zetac_eqs!(du, u, p)
    Δ, μ, ζ = u
    Δ = abs(Δ)
    β, int = p

    du[1] = A_int(β, μ, ζ, Δ, int)[1]
    du[2] = 1/(3π^2)-numberDensity_mf(β, μ, ζ, Δ)[1]
    du[3] = Omega_sp(β, μ, ζ, Δ, int)[1] - Omega_sp(β, μ, ζ, 0.0, int)[1]
    return nothing
 end

 function find_zetac_Deltac_muc(int, β = 1e3, zeta0 = 0.1, solve_alg = LevenbergMarquardt())
        Δ0, μ0 = find_Δ_μ(min(1e3, β), int, 0.0) # Quite sensitive to initial guess; imbalanced version is worse, for some reason.

        prob = NonlinearProblem(zetac_eqs!,[Δ0, μ0, zeta0],[β, int])
        sol = solve(prob,solve_alg)#, maxiters = 100000)

        if !SciMLBase.successful_retcode(sol.retcode)
            @warn "The solver did not converge. Retcode: $(sol.retcode)"
            return NaN, NaN, NaN
        end
        return (abs(sol[3]), abs(sol[1]), sol[2])
 end

find_zetac(int, β = 1e3, zeta0 = 0.1) = find_zetac_Deltac_muc.(int, β, zeta0)[1]

interactions = -2.0:0.05:3.0

critical_zeta1 = zeros(length(interactions))
chemical_potential0 = zeros(length(interactions))

guess = 0.03
 for j = 1:(length(interactions))
    int = interactions[j]

    if j ==1
        guess = 0.03
    else
        guess = critical_zeta1[j-1]
    end

    chemical_potential0[j] = find_Δ_μ(1e3, interactions[j])[2]
   
    critical_zeta1[j], Δ, μ = find_zetac_Deltac_muc(interactions[j], 1e3, guess)

    println("Critical zeta for int = $(interactions[j]): $(critical_zeta1[j]), guess = $(guess)")
    println("Δ = $(Δ)")
end

# the minimum of the bogoliubov dispersion gives the separation between the unpolarized phase and the polarized phase

min_bogoliubov_disp = zeros(length(interactions))


 for j = 1:(length(interactions))
    int = interactions[j]

    Δ, μ = find_Δ_μ(1e3, interactions[j], 0.0)
    
    min_bogoliubov_disp[j] = sqrt(Δ^2 + min(μ,0.0)^2)
end

 # Function to find separation between mixed phase and normal phase

 using NonlinearSolve, SciMLBase

 zeroifpositive(x) = x > 0 ? 0.0 : x
 
 function zetac2_eqs!(du, u, p)
    Δ, μ, ζ = u
    Δ = abs(Δ)
    β, int = p

    du[1] = A_int(β, μ, ζ, Δ, int)[1] + abs(zeroifpositive(G_int(β, μ, ζ, Δ)[1]))^0.5
    # Second term is added to make sure second derivative is positive
    du[2] = 1/(3π^2)-numberDensity_mf(β, μ, ζ, 0.0)[1]
    du[3] = Omega_sp(β, μ, ζ, Δ, int)[1] - Omega_sp(β, μ, ζ, 0.0, int)[1]
    return nothing
 end



 function find_zetac2_Deltac2_muc2(int, β = 1e3, zeta0 = 0.1, solve_alg = LevenbergMarquardt())
        Δ0, μ0 = find_Δ_μ(min(1e3, β), int, 0.0) # Quite sensitive to initial guess; imbalanced version is worse, for some reason.

        prob = NonlinearProblem(zetac2_eqs!,[Δ0, μ0, zeta0],[β, int])
        sol = solve(prob,solve_alg, reltol = 1e-4, abstol = 1e-6)#, maxiters = 100000)

        if !SciMLBase.successful_retcode(sol.retcode)
            @warn "The solver did not converge. Retcode: $(sol.retcode)"
            return NaN, NaN, NaN
        end
        return (abs(sol[3]), abs(sol[1]), sol[2])
 end

find_zetac2(int, β = 1e3, zeta0 = 0.1) = find_zetac2_Deltac2_muc2.(int, β, zeta0)[1]



interactions2 = -2.0:0.05:2.35#2.5

critical_zeta2 = zeros(length(interactions2))


guess = 0.03
 for j = 1:(length(interactions2))
    int = interactions2[j]

    if j ==1
        guess = 0.03
    elseif int == 0.1
        guess = 0.7
    elseif 0.75 <= int <= 0.95
        guess = 1.0
    elseif 1.7 <= int <= 1.75
        guess =4.5
    elseif int == 1.85
        guess = 5.0
    else
        guess = isnan(critical_zeta2[j-1]) ? 1.0 : critical_zeta2[j-1]
    end

    critical_zeta2[j], Δ, μ = find_zetac2_Deltac2_muc2(interactions2[j], 1e3, guess, NewtonRaphson())

    println("Critical zeta2 for int = $(interactions2[j]): $(critical_zeta2[j])")
end


critical_zeta_low = min.(critical_zeta1, min_bogoliubov_disp)

using DelimitedFiles

writedlm("$PROJECT_ROOT/JLTP2026/critical_ζ.dat", hcat(interactions, critical_zeta_low), '\t')

# Make figure appendix

include("$(PROJECT_ROOT)/JLTP2026/initialize_makie.jl")

f = Figure(size = (490pt,240pt), figure_padding = (0,0,0,1))
ax1 = Axis(f[1, 2], xlabel = L"(k_F a_s)^{-1}", ylabel = L"\Delta/E_F", xticks = -2:1:3, limits = (-2,3.2,0,2.35))

zetaColors = [ColorSchemes.tab10[j] for j in 1:length(zetaVals)]

for j = 1:length(zetaVals)
    lines!(ax1, intVals, DeltaVals[:,j], label = L"ζ/E_F = %$(round(zetaVals[j], digits=2))", color = zetaColors[j])
end

critical_int = zeros(length(zetaVals)-1)

Legend(f[1,3], ax1)


ax2 = Axis(f[1, 1], xlabel = L"(k_F a_s)^{-1}", ylabel = L"ζ/E_F", xticks = -2:1:3, yticks = 0:2:10, limits=(-2,3,0,10))

nanif(x, condition) = condition ? x : NaN


lines!(ax2, interactions, nanif.(critical_zeta1, 1.0 .<= interactions .<= 2.35))
lines!(ax2, interactions2, critical_zeta2) 

lines!(ax2, interactions, nanif.(critical_zeta1, 1.0 .>= interactions ))
lines!(ax2, interactions, nanif.(critical_zeta1, interactions .>=  2.35 ))

lines!(ax2, interactions, nanif.(min_bogoliubov_disp, interactions .>= 1.0))

scatter!(ax2, [2.368], [6.876], color = :red, markersize = 8) # label = "Tricritical point" (Parish 2007)


text!(-0.5, 7.5, text = L"\text{N}", align = (:center, :center), fontsize = 14)
text!(2.5, 2.5, text = L"\text{SF}_0", align = (:center, :center), fontsize = 14)
text!(2.2, 7.5, text = L"\text{SF}_P \rightarrow", align = (:center, :center), fontsize = 14)
text!(0.92, 2.5, text = L"\text{PS}\rightarrow", align = (:center, :center), fontsize = 14)

Label(f[1,1,TopLeft()], "(a)", halign = :left)
Label(f[1,2,TopLeft()], "(b)", halign = :left)


save("$(PROJECT_ROOT)/JLTP2026/figures/fig5_zeta_saddlepoint.pdf", f)