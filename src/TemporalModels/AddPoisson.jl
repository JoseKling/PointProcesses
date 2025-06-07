export AddPoisson

struct AddPoisson <: ParametricTemporalModel
    I::Interval
    f
    ∫f
    ∫If::Real

    AddPoisson(args...) = initialize(args...)
end

function simulate(model::AddPoisson, params::Parameters{2})
    return sort!([simulate(params[1], model.I);
        thinning(simulate(params[2], model.I), model.f)])
end

function estimate(model::AddPoisson, events::Times)
    isempty(events) && return (0.0, 0.0)
    f_times   = model.f.(events)
    initial_x = [length(f_times) / measure(model.I), length(f_times) / measure(model.I)]
    lower     = [0.0, 0.0]
    upper     = [Inf, Inf]
    results   = optimize(x -> objective_ap(x, f_times, model.∫If, measure(model.I)),
                         (stor, x) -> gradient_ap!(stor, x, f_times, model.∫If, measure(model.I)),
                         (stor, x) -> hessian_ap!(stor, x, f_times),
                         lower,
                         upper,
                         initial_x,
                         IPNewton())
    return (Optim.minimizer(results)[1], Optim.minimizer(results)[2])
end

function objective_ap(x, f_times, ∫f, M)
    return (x[1] .* M) + (x[2] .* ∫f) - # integral of the CIF 
           sum(log.(x[1] .+ (x[2] .* f_times))) # Sum of log of f at event times
end

function gradient_ap!(storage, x, f_times, ∫f, M)
    storage[1] = M - sum(1.0 ./ (x[1] .+ (x[2] .* f_times)))
    storage[2] = ∫f - sum(f_times ./ (x[1] .+ (x[2] .* f_times)))
end

function hessian_ap!(storage, x, f_times)
    storage[1, 1] = sum(1.0 ./ ((x[1] .+ (x[2] .* f_times))) .^ 2.0)
    storage[1, 2] = sum(f_times ./ ((x[1] .+ (x[2] .* f_times))) .^ 2.0)
    storage[2, 1] = storage[1, 2]
    storage[2, 2] = sum((f_times .^ 2) ./ ((x[1] .+ (x[2] .* f_times))) .^ 2.0)
end

function rescaling!(model::AddPoisson, params::Parameters{2}, events::Times)
    T           = rescaling!(params[1], events, model.I)
    @. events .+= (params[2] * model.∫f.(model.I[1], events))
    return T + (params[2] * model.∫If)
end

CIF(model::AddPoisson, params::Parameters{2}, _...) = x -> params[1] + (params[2] * model.f(x))

∫CIF(model::AddPoisson, params::Parameters{2}, _...) =
    I -> (params[1] * (I[2] - I[1])) + (params[2] * model.∫f(I))

function min_max_CIF(model::AddPoisson, params::Parameters{2}, _...)
    return params[1], params[1] + params[2]
end



# function simulate!(model::AddPoisson, params::Parameters{2}, sim::Times)
#     n_const = simulate!(params[1], model.I, sim)
#     n_f     = simulate!(params[2], model.I, (@view sim[n_const+1: end]))
#     n_f     = thinning!((@view sim[n_const+1: n_const+n_f]), model.f)
#     sort!((@view sim[1: n_const+n_f]))
#     return n_const + n_f
# end

