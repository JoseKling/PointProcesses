export simulate, estimate, rescaling
export CIF, ∫CIF

"""
Abstract type for parametric models.
Functions expected to be implemented:

- Minimum (enough for residual based goodness-of-fit tests)
    - 'simulate'
    - 'estimate'
    - 'CIF'
- For better performance in some tests
    - '∫CIF'
- For rescaling based goodness-of-fit tests for temporal processes
    - '∫CIF' or
    - 'rescaling' (for better performance)

If '∫CIF' is not implemented, approximations using 'CIF' will be used.
If 'rescaling' is not implemented, will use '∫CIF'.
"""
abstract type ParametricModel end

"""
    simulate(model::ParametricModel, params::Parameters)

Simulates events from the model given the parameters. Number of parameters should
match the model's requirements. The function returns a sorted array of events.
"""
function simulate(model::ParametricModel, params::Parameters{n}) where n
    error("Simulation not implemented for $(typeof(model))")
end

"""
    estimate(model::ParametricModel, events::Events)

Estimates the parameters of the model based on the observed events. The function
returns the estimated parameters as a tuple.
"""
function estimate(model::ParametricModel, events::Events)
    error("Estimation not implemented for $(typeof(model))")
end

function CIF(model::ParametricModel, params::Parameters{n}, events...) where n # events might not be necessary
    error("CIF not implemented for $(typeof(model))")
end

function ∫CIF(model::ParametricModel, params::Parameters{n}, events...) where n # events might not be necessary
    @warn "∫CIF not implemented. Using default implementation."
    cif = CIF(model, params, events...)
    return integral(cif)
end

"""
    rescaling!(model::ParametricModel, params::Parameters, events::Events)

In-place re-scaling \$t \\mapsto \\int_0^t \\lambda(s) ds\$, where \$\\lambda(s)\$ is the
conditional intensity function. The resulting process is a homogeneous Poisson process.
"""
function rescaling!(model::ParametricModel, params::Parameters{n}, events::Events) where n
    error("rescaling! not implemented for $(typeof(model))")
end

