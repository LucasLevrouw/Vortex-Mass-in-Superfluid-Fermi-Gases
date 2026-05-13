## SYMBOLS

# Δ is the absolute value of Ψ
# β = 1/(kB*T)
# int = 1/(kF*a_s)


## UNITS

# 2m = 1



## BASIC DEFINITIONS


xi(k, μ) = k^2-μ # 2m = 1
En(k, μ, Δ) = sqrt(xi(k,μ)^2+Δ^2)

# functions defined as in Marini et al 1998
x1(x0) = sqrt((sqrt(1+x0^2) + x0) / 2)
kappa(x0) = x1(x0)^2/(1+x0^2)^(1/2) # defined without the square -> convention elliptic functions
I5(x0) = (1+x0^2)^(1/4) * ellipe(kappa(x0)) - 1/(2x1(x0)^2) * I6(x0)
I6(x0) = 1/(2(1+x0^2)^(1/4))  * ellipk(kappa(x0))





function En_minus_xi(k::Real, μ::Real, Δ::Real)
   if k^2-μ ≤ 0
       return En(k,μ,Δ) - xi(k,μ)
   else
       return Δ^2 / (En(k,μ,Δ) +xi(k,μ))
   end
end



function gapEqIntegrand(k::Real, β::Real, μ::Real, ζ::Real, Δ::Real) # == 1 - 2k^2 * f1(En(k, μ, Δ), β, ζ)
    ϵ = En(k, μ, Δ)
    ζ = abs(ζ)

    if β == Inf && ϵ > ζ
        # 1 - k^2 * 1/ϵ == (ϵ - k^2) / ϵ
        retval =  (μ^2 + Δ^2 - 2*k^2*μ) / (ϵ * (ϵ + k^2))
    elseif ϵ < ζ
        return 1 - 2k^2 * f1(ϵ, β, ζ)
    else
        if ϵ ≥ 1e100
            retval = 1.0 - 2k^2 
        else
        retval =  1/(1 + exp(-2β*ϵ) +
            exp(-β*(ϵ+ζ))+ exp(-β*(ϵ-ζ)))*
            ((μ^2 + Δ^2 - 2*k^2*μ)/(ϵ*(ϵ + k^2)) +
            (1 + k^2/ϵ)*exp(-2β*ϵ) + exp(-β*(ϵ+ζ))+ exp(-β*(ϵ-ζ)))
        end
    end
    if isnan(retval) || isinf(retval)
        @warn "gapEqIntegrand is $retval at k = $k, β = $β, μ = $μ, ζ = $ζ, Δ = $Δ; ϵ ≥ ζ: $(ϵ ≥ ζ)"
    end

    return retval
end

function A_int(β::Real, μ::Real, ζ::Real, Δ::Real, int::Real)
    # println("Calculating A_int for β = $β, μ = $μ, ζ = $ζ, Δ = $Δ")
    function integrand_(k)
        integrand =  gapEqIntegrand(k, β, μ, ζ, Δ)
        if isnan(integrand)
            @warn "Integrand is NaN at k = $k, β = $β, μ = $μ, ζ = $ζ, Δ = $Δ"
        end

        return integrand
    end
    integral, error = quadgk(k -> integrand_(k), 0, Inf)
    # integral, error = quadgk(k -> gapEqIntegrand(k, β, μ, ζ, Δ), 0, Inf)
    return -int/(8*π) + 1/(4*π^2)*integral[1], 1/(4*π^2)*error[1]
end

function numberDensityIntegrand(k::Real, β::Real, μ::Real, ζ::Real, Δ::Real) # == k^2 * (1 - 2xi(k,μ) * f1(En(k, μ, Δ), β, ζ)) 
    ϵ = En(k, μ, Δ)
    ζ = abs(ζ)
    if ϵ >= 1e100
        return k^2 * (1 - 2xi(k,μ))
    end

    if β == Inf
        if ϵ ≥ ζ
            retval =  k^2 /ϵ * En_minus_xi(k,μ,Δ) # zero temperature limit
        else
            retval = k^2  # zero temperature limit
        end
    else
        if ϵ ≥ ζ
            retval = 1/(1 + exp(-2β*ϵ) + exp(-β*(ϵ-ζ))+exp(-β*(ϵ+ζ)))*(k^2 /ϵ * En_minus_xi(k,μ,Δ) +
            (ϵ + xi(k, μ))/ϵ * k^2 * exp(-2β*ϵ) + k^2*exp(-β*(ϵ-ζ)) + k^2*exp(-β*(ϵ+ζ)))
        else
            # 1/(1 + exp(-2β*ϵ) + exp(-β*ϵ)*(exp(β*ζ)+exp(-β*ζ)))*(k^2 /ϵ * En_minus_xi(k,μ,Δ) +
            # (ϵ + xi(k, μ))/ϵ * k^2 * exp(-2β*ϵ) + 2*k^2*exp(-β*ϵ)*cosh(β*ζ))
            # I could rewrite the above; but this might be good enough, since for k->infty we're in the other case anyway
            retval = k^2 * (1 - 2xi(k,μ) * f1(ϵ, β, ζ)) 
        end
    end
    if isnan(retval) || isinf(retval)
        @warn "numberDensityIntegrand is $retval at k = $k, β = $β, μ = $μ, ζ = $ζ, Δ = $Δ; ϵ ≥ ζ: $(ϵ ≥ ζ)"
    end
    return retval
end


# returns the density (and the error)
numberDensity_mf(β, μ, ζ, Δ) =
    1/(2π^2).*quadgk(k -> numberDensityIntegrand(k, β, μ, ζ, Δ), 0, Inf)

# numberEqRHS_mf(β, μ, ζ, Δ) = numberDensity_mf(β, μ, ζ, Δ)[1]

function imbalanceDensityIntegrand(k::Real, β::Real, μ::Real, ζ::Real, Δ::Real) # == k^2 * sinh(β*ζ) / (cosh(β*ζ) + cosh(β*En(k, μ, Δ)))
    return k^2 * X(ζ, β, En(k, μ, Δ))
    # if β == Inf
    #     return k^2 * (ζ > En(k, μ, Δ))
    # elseif ζ > ϵ
    #     return k^2 * (1 + exp(-2β*ζ)) / (exp(-β*ϵ - β*ζ) + exp(β*ϵ - β*ζ) + exp(-2β*ζ) + 1)
    # else
    #     return k^2 * (1 - exp(-2β*ϵ)) / (1 + exp(-2β*ϵ) + exp(-β*ϵ - β*ζ) + exp(-β*ϵ + β*ζ))
    # end
end
function imbalanceDensity_mf(β::Real, μ::Real, ζ::Real, Δ::Real)
    return 1/(2π^2) .* quadgk(k -> imbalanceDensityIntegrand(k, β, μ, ζ, Δ), 0, Inf)
end

    function Omega_sp(β::Real, μ::Real, ζ::Real, Δ::Real, int::Real)
   # integral, error = 1/(2π^2).*quadgk(k-> 
   #    k^2*((1/β)*log(1 + exp(-2β*En(k, μ, Δ)) +
   #    2exp(-β*En(k, μ, Δ))*cosh(β*ζ))) + 
   #    Δ^2/(2(En(k, μ, Δ) + xi(k, μ))) *
   #    ((2μ*k^2 - μ^2 - Δ^2)/(En(k, μ, Δ) + k^2) +
   #    (2k^2*μ - μ^2)/(xi(k, μ) + k^2)) , 0, Inf)

    if β == 0.0
        return Inf
    end

    # println("Calculating Omega_sp for β = $β, μ = $μ, ζ = $ζ, Δ = $Δ")
    integral, error = 1/(2π^2).*quadgk(k-> Omega_sp_integrand(k, β, μ, ζ, Δ)
     , 0, Inf)

   return -int*Δ^2/8π - integral, error
end



function Omega_sp_integrand(k::Real, β::Real, μ::Real, ζ::Real, Δ::Real)
# == k^2 * ( (1/β)*log(2cosh(β*En(k, μ, Δ)) + 2cosh(β*ζ)) - xi(k,μ) - Δ^2/(2k^2) )
    ζ = abs(ζ)
    ϵ = En(k, μ, Δ)
    if ϵ >= 1e100
        retval = k^2*En_minus_xi(k, μ, Δ) - Δ^2/2 
    elseif β == Inf
        if ϵ ≥ ζ
            retval = k^2*En_minus_xi(k, μ, Δ) - Δ^2/2 
        else
            retval = k^2*(ζ - xi(k,μ)) - Δ^2/2 
        end
    else
        if ϵ ≥ ζ
            if k^2-μ ≤ 0
                retval = k^2*((1/β)*log(1 + exp(-2β*ϵ) + exp(-β*ϵ+β*ζ) + exp(-β*ϵ-β*ζ))) + 
                k^2*(En(k, μ, Δ)-xi(k,μ)) - Δ^2/2 
            else
                retval = k^2*((1/β)*log(1 + exp(-2β*ϵ) + exp(-β*ϵ+β*ζ) + exp(-β*ϵ-β*ζ))) + 
                # k^2 * Δ^2 / (En(k,μ,Δ) +xi(k,μ)) - Δ^2/2 
                # -Δ^2/2 * (En(k, μ, Δ) + xi(k,μ) - 2k^2)/(En(k, μ, Δ)  + xi(k,μ))
                -Δ^2/2 * 1/(En(k, μ, Δ)  + xi(k,μ)) * (-μ + (Δ^2 + μ^2 - 2μ*k^2)/(En(k, μ, Δ) + k^2))
            end
        else
            retval = k^2*((1/β)*log(1 + exp(-2β*ζ) + exp(-β*ζ+β*ϵ) + exp(-β*ζ-β*ϵ))) + 
            k^2*(ζ - xi(k,μ)) - Δ^2/2 
        end
    end
    if isnan(retval) || isinf(retval)
        @warn "Omega_sp integrand is $retval at k = $k, β = $β, μ = $μ, ζ = $ζ, Δ = $Δ; ϵ ≥ ζ: $(ϵ ≥ ζ)"
    end
    return retval
end


using VortexMass: I5, I6
using StaticArrays
using NLsolve
using NonlinearSolve, SciMLBase


function μN_eq(du, u, p)
    μ = u[1]
    β, ζ = p
    du[1] = 1/(3π^2) - numberDensity_mf(β, μ, ζ, 0.0)[1]
end

function find_μ_normal_state(β::Real, ζ::Real)
    sol = nlsolve((du,u)-> μN_eq(du, u, (β,ζ)), [0.9])
    return sol.zero[1]
end

function saddle_point_eqs(u,p)
    Δ, μ = u
    β, ζ, int = p
    # if μ > 1 return [1000,1000] end
    return SA[A_int(β, μ, ζ, Δ, int)[1], 1/(3π^2)-numberDensity_mf(β, μ, ζ, Δ)[1]]
end

function saddle_point_eqs!(du, u,p)
    Δ, μ = u
    β, ζ, int = p
    # if μ > 1 return [1000,1000] end
    du[1] = A_int(β, μ, ζ, Δ, int)[1]
    du[2] = 1/(3π^2)-numberDensity_mf(β, μ, ζ, Δ)[1]
end

function find_Δ_μ(β::Real, int::Real, ζ::Real = 0.0, solver = NewtonRaphson())
    # Use these as starting values for the solver
    x0_T0 = 1/8*exp(-π/2*int+2) - sqrt( (int > 0) * 3π/16*int^3 ) # approximate value at T = 0
    Δ_T0 = 1/(x0_T0*I5(x0_T0)+I6(x0_T0))^(2/3)
    μ_T0 = x0_T0 * Δ_T0

    # solve saddle point equations
    # prob = NonlinearProblem(saddle_point_eqs,SA[Δ_T0, μ_T0],SA[β, ζ, int])
    prob = NonlinearProblem(saddle_point_eqs!,[Δ_T0, μ_T0],[β, ζ, int])
    sol = solve(prob,solver)

    # if !SciMLBase.successful_retcode(sol.retcode)
    #     @warn "The solver did not converge. Returning normal state solution."
    #        μN = find_μ_normal_state(β, ζ)
    #     return (0.0, μN)
    # else
    #     return (abs(sol[1]), sol[2])
    # end

    if Omega_sp(β, sol[2], ζ, 0.0, int)[1] < Omega_sp(β, sol[2], ζ, sol[1], int)[1]
        μN = find_μ_normal_state(β, ζ)
        return (0.0, μN)
    else    
        if !SciMLBase.successful_retcode(sol.retcode)
            @warn "The solver did not converge."
        end
        return (sol[1], sol[2])
    end
end


function Tc_eqs!(du, u,p)
   x, μ = u
   int, Δ = p
   ζ = 0.0
   # if μ > 1 return [1000,1000] end
   du[1] = A_int(2^x, μ, ζ, 0.0, int)[1]
   du[2] = 1/(3π^2)-numberDensity_mf(2^x, μ, ζ, 0.0)[1]
end

function find_Tc(int::Real, equation_of_state="MF")
   find_Tc_μc(int, equation_of_state)[1]
end

function find_Tc_μc(int::Real, equation_of_state="MF")
   # β = 2^x
   if equation_of_state == "MF" || equation_of_state == "SP"
       # prob = NonlinearProblem(Tc_eqs!,[1.0, 0.5],int)
       # sol = solve(prob, NewtonRaphson(autodiff=false))
       sol = nlsolve((du,u)-> Tc_eqs!(du,u,(int,0.0)), [1.0, 0.9])
       if !sol.x_converged && !sol.f_converged
           @warn "Not converged."
       end 
       # print error on equations with found solution 
      #  res = [0.0,0.0]
      #  Tc_eqs!(res, sol.zero,(int,0.0))
      #  println(res)
       return [2^(-sol.zero[1]), sol.zero[2]]
   elseif equation_of_state == "GPF"
           throw("not implemented yet")
   else
       throw(ArgumentError("""provide a valid value for `equation_of_state`"""))
   end
end

function X(ϵ::Real, β::Real, ζ::Real) # == sinh(β*ϵ) / (cosh(β*ζ) + cosh(β*ϵ))
     if !(ϵ ≥ 0  && β ≥ 0)
        @warn "X is not defined for ϵ = $ϵ, β = $β, ζ = $ζ. Returning NaN."
    end

    @assert ϵ ≥ 0  && β ≥ 0
   
     ζ = abs(ζ)
    if β == Inf
        return ϵ > abs(ζ) ? 1.0 : 0.0
    else
        return ϵ > abs(ζ) ? (1 - exp(-2β*ϵ))/ (1 + exp(-2β*ϵ) + exp(-β*ϵ+β*ζ)+exp(-β*ϵ-β*ζ)) : (exp(β*(ϵ-abs(ζ))) - exp(-β*(ϵ+abs(ζ))))/ (1 + exp(-2β*abs(ζ)) + exp(β*(ϵ-abs(ζ)))+ exp(-β*(ϵ+abs(ζ))))
    end
end

function Y(ϵ::Real, β::Real, ζ::Real) # == (∂/∂ϵ) X(ϵ, β, ζ) == β (1 + cosh(β*ϵ) cosh(β*ζ)) / (cosh(β*ϵ) + cosh(β*ζ))^2
    ζ = abs(ζ)
    @assert ϵ ≥ 0  && β ≥ 0 
    if β == Inf
        error("Y is not defined at zero temperature")
    elseif ϵ < ζ
        return Y(ζ, β, ϵ) # Y is symmetric in ϵ and ζ
    else
    # return 2β * (2exp(-2β*ϵ) +cosh(β*ζ)*( exp(-β*ϵ) +exp(-3β*ϵ)))/ (1 + exp(-2β*ϵ) + 2exp(-β*ϵ)*cosh(β*ζ))^2
        return β * (4exp(-2β*ϵ) + exp(-β*(ϵ+ζ)) + exp(-β*(ϵ-ζ)) +exp(-3β*ϵ+β*ζ) + exp(-3β*ϵ-β*ζ)) /
                (1 + exp(-2β*ϵ) + exp(-β*(ϵ+ζ)) + exp(-β*(ϵ-ζ)) )^2
    end
end

f1(ϵ, β, ζ) =  X(ϵ, β, ζ) / (2ϵ)

f2(ϵ, β, ζ) = 1/(4ϵ^2) * (X(ϵ, β, ζ)/ϵ - Y(ϵ, β, ζ)) # == -1/(2ϵ) (∂/∂ϵ) f1(ϵ, β, ζ) 

function f3(ϵ::Real, β::Real, ζ::Real)  # == -1/(4ϵ) (∂/∂ϵ) f2(ϵ, β, ζ)
    ζ = abs(ζ) 
    @assert ϵ ≥ 0  && β ≥ 0 
    if β == Inf
        @warn "f3 is not defined at zero temperature"
        return NaN
    elseif ζ > ϵ
        f3num1 = 1/4 * (-3β*ϵ * (1 + exp(-4β*ζ) + 6exp(-2β*ζ)) * (exp(-β*ϵ - β*ζ) + exp(β*ϵ - β*ζ))  + 
        2(-exp(-β*ϵ - β*ζ) + exp(β*ϵ - β*ζ)) * (3/2 * (exp(-2β*ϵ - 2β*ζ) + exp(2β*ϵ - 2β*ζ)) + 
        exp(-2β*ζ) * (6 - 3β^2*ϵ^2) + 1/2*(1 + exp(-4β*ζ)) * (3 +β^2*ϵ^2)) -
        (1 + exp(-2β*ζ)) * (18β*ϵ * exp(-2β*ζ) + exp(-2β*ϵ - 2β*ζ) * (6 + 3β*ϵ -β^2*ϵ^2) + 
        exp(2β*ϵ - 2β*ζ) * (-6 + 3β*ϵ +β^2*ϵ^2)))

        f3den1 = 4ϵ^5 * (exp(-β*ϵ - β*ζ) + exp(β*ϵ - β*ζ) + exp(-2β*ζ) + 1)^3 
        return f3num1 / f3den1
    else
        f3num2= 1/4 * (-3β*ϵ * (exp(-2β*ϵ) + 1) * (6exp(-2β*ϵ) + exp(-2β*ϵ - 2β*ζ) + exp(-2β*ϵ + 2β*ζ)) + 
        2(-exp(-2β*ϵ) + 1) * ((6 - 3β^2*ϵ^2)*exp(-2β*ϵ) + 3/2 * (exp(-4β*ϵ) + 1) + 1/2 * (exp(-2β*ϵ - 2β*ζ) + 
        exp(-2β*ϵ + 2β*ζ)) * (3+ β^2*ϵ^2)) - (exp(-β*(ϵ + ζ)) + exp(-β*(ϵ - ζ))) * ((6 + 3β*ϵ - β^2*ϵ^2)*exp(-4β*ϵ) + 
        18β*ϵ * exp(-2β*ϵ) + (-6 + 3β*ϵ + β^2*ϵ^2)))

        f3den2 = 4ϵ^5 * (1 + exp(-2β*ϵ) + exp(-β*ϵ - β*ζ) + exp(-β*ϵ + β*ζ))^3
        return f3num2 / f3den2
    end
end 

function f4(ϵ::Real, β::Real, ζ::Real) # == -1/(6ϵ) (∂/∂ϵ) f3(ϵ, β, ζ) 
    ζ = abs(ζ)
    @assert ϵ > 0 && β ≥ 0 
    if β == Inf
        error("f4 is not defined at zero temperature")
    elseif ζ > ϵ
        f4num1 = 
        1/2* (3(-exp(-β*(ϵ + ζ)) + exp(-β*(-ϵ + ζ))) * (1 + exp(-6β*ζ)) * (5 + 2β^2*ϵ^2) - 
        β*ϵ*(exp(-β*ϵ - β*ζ) + exp(β*ϵ - β*ζ)) * (exp(-2β*ζ) + 1) *
        ((195 - 11β^2*ϵ^2)*exp(-2β*ζ) + (15 +  β^2*ϵ^2) * (exp(-4β*ζ) + 1) )
        + (1 + exp(-4β*ζ)) * (45(-exp(-2β*(ϵ + ζ)) + exp(-2β*(-ϵ + ζ))) - 
        8β*ϵ*(15 + β^2*ϵ^2)*exp(-2β*ζ)) + 2β*ϵ*(exp(-2β*ϵ - 2β*ζ) + exp(2β*ϵ - 2β*ζ)) *
        (-4(15 +   β^2*ϵ^2)*exp(-2β*ζ) + (exp(-4β*ζ) +  1) * (-15 + 2β^2*ϵ^2)) + 
        6(5/2*  (-exp(-4β*ϵ - 4β*ζ) + exp( 4β*ϵ - 4β*ζ)) + (-60β*ϵ +   8β^3*ϵ^3)*exp(-4β*ζ) - 
        4(-5 + 2β^2*ϵ^2) * (-exp(-2β*ϵ - 4β*ζ) + exp(2β*ϵ - 4β*ζ)) ) -
        (1 + exp(-2β*ζ)) * ((45 + 15β*ϵ - 6β^2*ϵ^2+ β^3*ϵ^3)*exp(-3β*(ϵ + ζ)) + 
        exp(-3β*(-ϵ + ζ)) * (-45 + 15β*ϵ +  6β^2*ϵ^2+ β^3*ϵ^3) + 
        18(5 - 2β^2*ϵ^2)*exp(-3β*ζ - β*ϵ) + 18(-5 + 2β^2*ϵ^2)*exp(-3β*ζ + β*ϵ)))

        f4den1 = 48ϵ^7 *(1 + exp(-2β*ζ) + exp(-β*ϵ - β*ζ) + exp(β*ϵ - β*ζ))^4
        return f4num1 / f4den1
    else
        f4num2=
        1/2*  (3(1 - exp(-2β*ϵ)) * (exp(-3β*(ϵ + ζ)) + exp(-3β*(ϵ - ζ))) * (5 + 2β^2*ϵ^2) - 
        β*ϵ*(exp(-2β*ϵ) + 1) * (exp(-β*ζ - β*ϵ) + exp(β*ζ - β*ϵ)) *
        ((195 - 11β^2*ϵ^2)*exp(-2β*ϵ) +(15 + β^2*ϵ^2) * (exp(-2β*ϵ - 2β*ζ) + exp(-2β*ϵ + 2β*ζ)) ) +
        (exp(-2β*(ϵ + ζ)) + exp(-2β*(ϵ - ζ))) * (45(1 - exp(-4β*ϵ)) - 
        8β*ϵ*(15+ β^2*ϵ^2)*exp(-2β*ϵ)) +  4β*ϵ*(exp(-4β*ϵ) + 1) * (-2(15+ β^2*ϵ^2)*exp(-2β*ϵ) + 
        1/2* (-15+ 2β^2*ϵ^2) * (exp(-2β*ϵ - 2β*ζ) + exp(-2β*ϵ + 2β*ζ)) ) + 
        6(5/2*  (-exp(-8β*ϵ) + 1) + (-60β*ϵ + 8β^3*ϵ^3)*exp(-4β*ϵ) - 
        4(-5 + 2β^2*ϵ^2)*exp(-2β*ϵ) * (-exp(-4β*ϵ) + 1) ) - (exp(-β*(ϵ + ζ)) +  exp(-β*(ϵ - ζ))) *
        ((45 + 15β*ϵ - 6β^2*ϵ^2+ β^3*ϵ^3)*exp(-6β*ϵ) + (-45 + 15β*ϵ + 6β^2*ϵ^2+ β^3*ϵ^3) + 
        18(-5 + 2β^2*ϵ^2) * (exp(-2β*ϵ) - exp(-4β*ϵ))))

        f4den2 = 48ϵ^7*(1 + exp(-2β*ϵ) + exp(-β*ϵ + β*ζ) + exp(-β*ϵ - β*ζ))^4
        return f4num2 / f4den2
    end
end 


Dt_int(β, μ, ζ, Δ) = quadgk(k-> 1/(2π^2) * k^2 * xi(k,μ) * f2(En(k,μ,Δ),β,ζ),
    0, Inf)


C_int(β, μ, ζ, Δ) = 1/(3π^2).*quadgk(k -> k^4*f2(En(k, μ, Δ), β, ζ), 0, Inf)

# function C_int(β, μ, ζ, Δ)
#     function integrand_(k)
#         if isnan(k^4 * f2(En(k, μ, Δ), β, ζ))
#             println("oh no! k = $k, β = $β, μ = $μ, ζ = $ζ, Δ = $Δ, En = $(En(k, μ, Δ)), f2 = $(f2(En(k, μ, Δ), β, ζ))")
#         end
#         return k^4 * f2(En(k, μ, Δ), β, ζ)
#     end
#     # println("Calculating C_int for β = $β, μ = $μ, ζ = $ζ, Δ = $Δ")
#     integral, error = quadgk(k -> integrand_(k), 0, Inf)
#     return 1/(3π^2)*integral, 1/(3π^2)*error
# end

E_int(β, μ, ζ, Δ) = 2/(3*π^2) .* quadgk( k -> k^4*xi(k, μ)^2*f4(En(k, μ, Δ), β, ζ), 0, Inf)
Et_int(β, μ, ζ, Δ) = 2 .*E_int(β, μ, ζ, Δ)

Q_int(β, μ, ζ, Δ) = 1/(4*π^2*Δ^2) .*
  quadgk(k ->
   k^2*(f1(En(k, μ, Δ), β, ζ) - (En(k, μ, Δ)^2 + xi(k, μ)^2)*
   f2(En(k, μ, Δ),β, ζ)), 0, Inf)

R_int(β, μ, ζ, Δ) = (1/(4*π^2*Δ^2)) .*
   quadgk(k->
   k^2*((f1(En(k, μ, Δ), β, ζ) + (En(k, μ, Δ)^2 - 3*xi(k, μ)^2)*
   f2(En(k, μ, Δ), β, ζ))/(3*Δ^2) + 
      4/3*(xi(k, μ)^2 - 2*En(k, μ, Δ)^2)*f3(En(k, μ, Δ), β, ζ) + 
      2*En(k, μ, Δ)^2*Δ^2*f4(En(k, μ, Δ), β, ζ))
   ,  0,Inf)
Rt_int(β, μ, ζ, Δ) = 2 .*R_int(β, μ, ζ, Δ)


G_int(β, μ, ζ, Δ) = 1/(2*π^2) .* quadgk((k -> k^2 * f2(En(k, μ, Δ),β, ζ)), 0, Inf)

using BSplineKit


function calculate_mf_eft_parameters(β::Real, int::Real, ζ::Real = 0.0; max = 4, step = 0.002, calculate_interpolations = true, warn = true)
    Δ, μ = find_Δ_μ(β, int, ζ)

    if Δ == 0.0
        warn && @warn "Δ == 0. Cannot calculate EFT coefficients in normal state."
        return (; Δ, μ, A_itp=w->NaN, Dt_itp=w->NaN, Ct=NaN, Et=NaN, Q=NaN, Rt=NaN, G=NaN, ξ2=NaN)
    end

    Ct = C_int(β,μ,ζ,Δ)[1] # C == Ct if 2m == 1
    Q = Q_int(β,μ,ζ,Δ)[1]
    G = G_int(β,μ,ζ,Δ)[1]

    Et = Et_int(β,μ,ζ,Δ)[1]
    Rt = Rt_int(β,μ,ζ,Δ)[1]
    
    ξ2 = sqrt(2Ct /(G*Δ^2))

    if calculate_interpolations

    # Interpolation of A and Dt using BSPlineKit.jl
    Δsq_vals = Δ == 0 ? range(0., 0.1, step = step) : Δ^2 .*range(0.,max, step=step)
    basis = BSplineBasis(BSplineOrder(3), Δsq_vals);
    A_itp = approximate(Δsq -> A_int(β,μ,ζ,sqrt(Δsq),int)[1], basis)
    Dt_itp = approximate(Δsq -> Dt_int(β,μ,ζ,sqrt(Δsq))[1], basis)
    Ω_itp = approximate(Δsq -> Omega_sp(β,μ,ζ,sqrt(Δsq),int)[1], basis)


    return (; Δ, μ, Ct, Et, Q, Rt, G, A_itp, Dt_itp, Ω_itp, ξ2)

else

 return (; Δ, μ, Ct, Et, Q, Rt, G, ξ2)
end
end