abstract type BoundaryCondition end

# --- Dirichlet Boundary Condition ---

struct DirichletBC{T <: AbstractFloat} <: BoundaryCondition
    boundary_value::T
end

# Catches Int, Rational, and Irrationals (like pi), converting to Float64; defaults to 1.0
DirichletBC(val::Real=1.0) = DirichletBC(Float64(val))


# --- Neumann Boundary Condition ---

struct NeumannBC{T <: AbstractFloat} <: BoundaryCondition
    boundary_derivative::T
end

# Catches Int, Rational, and Irrationals (like pi), converting to Float64; defaults to 0.0
NeumannBC(val::Real=0.0) = NeumannBC(Float64(val))


# --- Helper Function for Parsing ---

function parse_bc(spec)
    # Handle standalone symbols (defaults to homogeneous/zero conditions)
    if spec === :dirichlet
        return DirichletBC()
    elseif spec === :neumann
        return NeumannBC()
        
    # Handle tuples with specific values
    elseif spec isa Tuple && length(spec) == 2
        bc_type, val = spec
        if bc_type === :dirichlet
            return DirichletBC(val)
        elseif bc_type === :neumann
            return NeumannBC(val)
        end
    end

    throw(ArgumentError("Unsupported specification: $spec. Use :dirichlet, :neumann, or a tuple like (:dirichlet, 1.5)."))
end
