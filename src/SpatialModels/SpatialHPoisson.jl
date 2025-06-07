export SpatialHPoisson

struct SpatialHPoisson{d} <: ParametricSpatialModel{d}
    A::Area{d}
end
SpatialHPoisson(args...; d::Int=2) = SpatialHPoisson(convert_interval(args..., d=d))

simulate(model::SpatialHPoisson, params::Parameters{1}) = simulate(params[1], model.A)

estimate(model::SpatialHPoisson{d}, events::Locations{d}) where d = (length(events) / measure(model.A),)

CIF(model::SpatialHPoisson, params::Parameters{1}, _...) = x -> params[1]

∫CIF(model::SpatialHPoisson, params::Parameters{1}, _...) = A -> params[1] * measure(A)

