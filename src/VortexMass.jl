# units: ħ == 2*m == 1
# EFT eom: im dΨ/dt = (-d^2/dx^2 + V(x) + g*|Ψ|^2) * Ψ


module VortexMass

greet() = print("Hello Physicist!")


## Calculation
using LinearAlgebra
using OrdinaryDiffEq
using BSplineKit # interpolation with b-splines
using Optim

# ## Visualisation
# using Plots
# using LaTeXStrings
# using Printf # for formatting strings
# using CairoMakie # alternative plotting package

using SpecialFunctions
using QuadGK
using NonlinearSolve


# include("./eft_coefficients.jl")
# include("./saddlepoint.jl")
include("./eft_coefficients_and_mean_field.jl")
export C_int, D_int, Q_int, R_int, A_int, Dt_int, Et_int, Rt_int, G_int, Omega_sp
export numberDensity_mf, imbalanceDensity_mf
export find_Δ_μ, calculate_mf_eft_parameters, find_Tc

include("boundary_conditions.jl")
export BoundaryCondition, DirichletBC, NeumannBC

include("./vortex-profile.jl")
export calculate_profile

include("./energy.jl")
export free_energy, free_energy_callback


# ## General functions

# include("./helper_functions.jl")
# export to_real, to_complex

# ## Functions for solving the differential equation

# include("./eom_2d_real.jl")
# export solve_eft_2d_real
# # export f_eft_2d!, f_eft_2d_real!

# include("./eom_2d_real_pin.jl")
# export solve_eft_2d_real_pin, solve_eft_2d_real_pin2

# include("./eom_2d_complex.jl")
# export solve_eft_2d_complex

# ## Steady state

# include("./steady.jl")
# export steady_state, interpolated_profiles

# ## Initial conditions

# include("./initial_condition.jl")
# export healing_length, coherence_length
# export u0_vortex, u0_two_vortices


end # module VortexMass


