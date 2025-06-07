export HPoisson

struct HPoisson <: ParametricTemporalModel
    I::Interval

    HPoisson(arg...) = new(convert_interval(arg...))
end

simulate(model::HPoisson, params::Parameters{1}) = simulate(params[1], model.I)

estimate(model::HPoisson, events::Times) = (length(events) / measure(model.I),)

rescaling!(model::HPoisson, param::Parameters{1}, events::Times) =
    rescaling!(param[1], events, model.I)

CIF(model::HPoisson, params::Parameters{1}, _...) = x -> params[1]

∫CIF(model::HPoisson, params::Parameters{1}, _...) = I -> params[1] * (I[2] - I[1])

min_max_CIF(model::HPoisson, params::Parameters{1}, _...) = params[1], params[1]



# simulate!(model::HPoisson, params::Parameters{1}, sim::Times) = 
#     simulate!(params[1], model.I, sim)
