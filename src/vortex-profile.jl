

# dirichlet on one side, neumann on the other
function f_profile!(du, u, p, t)
    # u = f = |Ψ(r)| / Δ  (normalized, u → 1 at infinity)

    (; R, N, l, Δ, A_itp, Ct, Et, bc) = p
    dr = R/N

    @inbounds begin
        # i = 1: u[1] = 0 fixed (vortex core, Dirichlet at r = 0)
        du[1] = 0.0
    end


    @inbounds for i in 2:(N-1)
        r = (i-1)*dr

        Δf = (u[i+1] + u[i-1] - 2u[i])/dr^2 + 1/r * (u[i+1] - u[i-1])/2dr
        Δfsq = (u[i+1]^2 + u[i-1]^2 - 2u[i]^2)/dr^2 + 1/r * (u[i+1]^2 - u[i-1]^2)/2dr

        A = A_itp(Δ^2 * u[i]^2)
        du[i] = Ct*(Δf-l^2/r^2*u[i]) - (A + Et*Δ^2*Δfsq)*u[i]
    end

    @inbounds begin
    i = N
    r = (N-1)*dr
    if bc isa DirichletBC
        u_R = bc.boundary_value
    elseif bc isa NeumannBC
        u_R = u[i] + dr * bc.boundary_derivative
    else
        error("Unsupported boundary condition: $(typeof(bc))")
    end

    Δf   = (u_R + u[i-1] - 2u[i])/dr^2 + 1/r * (u_R - u[i-1])/2dr
    Δfsq = (u_R^2 + u[i-1]^2 - 2u[i]^2)/dr^2 + 1/r * (u_R^2 - u[i-1]^2)/2dr

    A = A_itp(Δ^2 * u[i]^2)
    du[i] = Ct*(Δf-l^2/r^2*u[i]) - (A + Et*Δ^2*Δfsq)*u[i]
    end

    return nothing
end

# f_profile!(du, u, p) = f_profile!(du, u, p, 0.0)

# 8 min 20, int = -2 takes a long time
# neuman: 15 min (desktop)
# only 0.95 (thinkpad): 11 min

# neumann, R = 200, N = 200: 43 min (thinkpad)


function calculate_profile(β::Real,int::Real,ζ::Real, R::Real, N::Integer, t::Real=1e4, bc=:neumann, l::Integer=1; callback=nothing, saving_cache=NaN, solver=nothing, solve_kwargs=(; save_everystep=false))

   mf_eft_params = calculate_mf_eft_parameters(β,int,ζ)
   calculate_profile(mf_eft_params, R, N, t, bc, l; callback, saving_cache, solver, solve_kwargs)
end


function calculate_profile(mf_eft_params::NamedTuple, R::Real, N::Integer,  t::Real=1e4, bc=:neumann, l::Integer=1; callback=nothing, saving_cache=NaN, solver = nothing, solve_kwargs = (; save_everystep=false))

    (; Δ, A_itp, Ct, Et) = mf_eft_params

    bc = bc isa BoundaryCondition ? bc : parse_bc(bc)

    r = (R/N) .* collect(0:(N-1))

    p = (; R, N, l, Δ, A_itp, Ct, Et, bc, saving_cache)
    prob_profile_ode = ODEProblem(f_profile!, @.(tanh(r)), (0.0, t), p)

    if !isnothing(callback)
        solve_kwargs = (; solve_kwargs..., callback=callback)
    end

    sol_ode = solve(prob_profile_ode, solver; solve_kwargs...)

    return [r sol_ode.u[end]]
end


