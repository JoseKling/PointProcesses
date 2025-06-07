export SpatialAddPoisson

struct SpatialAddPoisson{d} <: ParametricSpatialModel{d}
    A::Area{d}
    f
    ∫f
    ∫Af::Float64
end
SpatialAddPoisson(args...; d::Int=2) = SpatialAddPoisson(initialize(args...; d=d)...)

function simulate(model::SpatialAddPoisson, params::Parameters{2})
    return [simulate(params[1], model.A);
        thinning(simulate(params[2], model.A), model.f)]
end

function estimate(model::SpatialAddPoisson{d}, events::Locations{d}) where d
    isempty(events) && return (0.0, 0.0)
    f_times   = model.f.(events)
    initial_x = [length(f_times) / measure(model.A), length(f_times) / measure(model.A)]
    lower     = [0.0, 0.0]
    upper     = [Inf, Inf]
    results   = optimize(x -> objective_sap(x, f_times, model.∫Af, measure(model.A)),
                         (stor, x) -> gradient_sap!(stor, x, f_times, model.∫If, measure(model.A)),
                         (stor, x) -> hessian_sap!(stor, x, f_times),
                         lower,
                         upper,
                         initial_x,
                         IPNewton())
    return (Optim.minimizer(results)[1], Optim.minimizer(results)[2])
end

function objective_sap(x, f_times, ∫f, M)
    return (x[1] .* M) + (x[2] .* ∫f) - # integral of the CIF 
           sum(log.(x[1] .+ (x[2] .* f_times))) # Sum of log of f at event times
end

function gradient_sap!(storage, x, f_times, ∫f, M)
    storage[1] = M - sum(1.0 ./ (x[1] .+ (x[2] .* f_times)))
    storage[2] = ∫f - sum(f_times ./ (x[1] .+ (x[2] .* f_times)))
end

function hessian_sap!(storage, x, f_times)
    storage[1, 1] = sum(1.0 ./ ((x[1] .+ (x[2] .* f_times))) .^ 2.0)
    storage[1, 2] = sum(f_times ./ ((x[1] .+ (x[2] .* f_times))) .^ 2.0)
    storage[2, 1] = storage[1, 2]
    storage[2, 2] = sum((f_times .^ 2) ./ ((x[1] .+ (x[2] .* f_times))) .^ 2.0)
end

CIF(model::SpatialAddPoisson, params::Parameters{2}, _...) = x -> params[1] + (params[2] * model.f(x))

∫CIF(model::SpatialAddPoisson, params::Parameters{2}, _...) =
    A -> (params[1] * prod(A)) + (params[2] * model.∫f(A))
