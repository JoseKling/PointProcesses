export Renewal

struct Renewal{D<:UnivariateDistribution} <: ParametricTemporalModel
    distribution::D
    function Renewal(dist::UnivariateDistribution)
        @assert minimum(dist) >= 0 "Distribution must have non-negative support"
        @assert maximum(dist) > 0 "Distribution must have positive support"
        new{typeof(dist)}(dist)
    end
end

function simulate(model::Renewal, region::Interval)
    times = region[1] .+ cumsum(rand(model.distribution, round(Int, 1.2 * measure(region) / mean(model.distribution))))
    times[end] > region[2] && return TPP(times[1:searchsortedfirst(times, region[2]) - 1], region)
    while true
        new_t = rand(model.distribution) + times[end]
        new_t > region[2] && break
        push!(times, new_t)
    end
    return TPP(times, region)
end

function estimate(model::Renewal, tpp::TPP)
    X = interevents(tpp)
    return Renewal(fit(typeof(model.distribution), X))
end

function intensity(model::Renewal, tpp::TPP)
    function cif(t)
        idx = searchsortedfirst(tpp.times, t) - 1
        ti = idx == 0 ? 0 : tpp.times[idx]
        return pdf(model.distribution, t - ti) / ccdf(model.distribution, t - ti)
    end
    return DomainFunction(cif, tpp.region)
end