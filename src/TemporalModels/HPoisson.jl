export HPoisson

struct HPoisson <: ParametricModel
    I::Interval

    HPoisson(I::Tuple{Real, Real}) = new(convert_interval(I))
    HPoisson(a::Real, b::Real) = HPoisson((a, b))
end

simulate(model::HPoisson, params::Parameters{1}) = simulate(params[1], model.I)

estimate(model::HPoisson, events::Events) = (length(events) / measure(model.I),)

time_transform(model::HPoisson, param::Parameters{1}, events::Events) =
    time_transform(param[1], events, model.I)

time_transform!(model::HPoisson, param::Parameters{1}, events::Events) =
    time_transform!(param[1], events, model.I)

CIF(model::HPoisson, params::Parameters{1}, _...) = x -> params[1]

∫CIF(model::HPoisson, params::Parameters{1}, _...) = I -> params[1] * (I[2] - I[1])

min_max_CIF(model::HPoisson, params::Parameters{1}, _...) = params[1], params[1]



# simulate!(model::HPoisson, params::Parameters{1}, sim::Events) = 
#     simulate!(params[1], model.I, sim)
