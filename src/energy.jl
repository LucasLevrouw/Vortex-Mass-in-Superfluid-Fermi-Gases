function free_energy(mf_eft_params, u, R::Real, bc = :neumann, l::Integer=1; use_trapz = true)

    (; Δ, Ct, Et, Ω_itp) = mf_eft_params

    E = Et/2

    bc = bc isa BoundaryCondition ? bc : parse_bc(bc)

    N = length(u)
    dr = R/N

    Ω_inf = Ω_itp(Δ^2)

    F = 0.0

    @inbounds for i in 2:(N-1)
        r = (i-1)*dr

        df_dr    = (u[i+1] - u[i-1])/(2dr)
        dfsq_dr  = (u[i+1]^2 - u[i-1]^2)/(2dr)

        F += 2π*r * (Ω_itp(Δ^2*u[i]^2) - Ω_inf + Ct*Δ^2 * (df_dr^2 + l^2/r^2*u[i]^2) - E*Δ^4 * dfsq_dr^2)
    end

    @inbounds begin
        # i = 1: r = 0, so the integrand vanishes — no contribution needed

        i = N
        r = (N-1)*dr

        u_R = bc isa DirichletBC ? bc.boundary_value : u[i] + dr * bc.boundary_derivative

        df_dr   = (u_R    - u[i-1])/(2dr)
        dfsq_dr = (u_R^2  - u[i-1]^2)/(2dr)

        FN = 2π*r * (Ω_itp(Δ^2*u[i]^2) - Ω_inf + Ct*Δ^2 * (df_dr^2 + l^2/r^2*u[i]^2) - E*Δ^4 * dfsq_dr^2)
        F += use_trapz ? FN/2 : FN
    end
    
    return F * dr
end

using DiffEqCallbacks

function affect!(integrator, mf_eft_params)
    R = integrator.p.R
    l = integrator.p.l
    bc = integrator.p.bc
    F = free_energy(mf_eft_params, integrator.u, R, bc, l)
    integrator.p = (integrator.p..., saving_cache = F)
end

function free_energy_callback(mf_eft_params, l::Integer=1)
    cb_update_free_energy = DiscreteCallback((u,t,integrator)-> true, integrator -> affect!(integrator, mf_eft_params); save_positions = (false, false))
    saved_values = SavedValues(Float64, Float64)
    cb_save_free_energy = SavingCallback((u, t, integrator) -> integrator.p.saving_cache, saved_values)

    cb=CallbackSet(cb_update_free_energy,cb_save_free_energy)
    return cb, saved_values
end