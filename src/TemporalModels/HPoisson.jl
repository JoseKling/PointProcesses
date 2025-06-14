simulate(model::HPoisson, region::Interval) = TPP(simulate(model.μ, region), region)

estimate(::HPoisson, tpp::TPP) = HPoisson(length(tpp.times) / measure(tpp.region))

rescale(model::HPoisson, tpp::TPP) = TPP(model.μ .* tpp.times, (0.0, measure(tpp.region) * model.μ))

intensity(model::HPoisson, tpp::TPP) = DomainFunction(x -> model.μ, tpp.region, integral= x -> model.μ * (x[2] - x[1]))

function simulate!(model::HPoisson, region::Interval, sim::Vector{Float64})
    return simulate!(model.μ, region, sim)
end

function rescale!(model::HPoisson, tpp::TPP, rescaled::Vector{Float64})
    @. rescaled[1:length(tpp.times)] = model.μ * tpp.times
end
