#---------------- Kolmogorov-Smirnov distance ---------------------------------
function statistic(model::ParametricTemporalModel, ::Type{KSExponential}, params::Parameters, events::Times)
    (length(events) <= 2) && return 1.0 # If 'events' has only 2 elements, than there is only one event, so no interevent times
    xs   = diff(rescaling!(model, params, events)[1])
    vals = ecdf(xs)
    d    = 0.0
    @inbounds for i in eachindex(xs)
        exp_cdf = 1.0 - exp.(-xs[i]) # The empirical CDF of the exponential distribution
        d       = max(d, exp_cdf - vals[i], vals[i+1] - exp_cdf)
    end
    return d
end

#--------------- Laplace distance ---------------------------------------------
function statistic(model::ParametricTemporalModel, ::Type{LPExponential}, params::Parameters, events::Times)
    (length(events) <= 2) && return 1.0 # If 'events' has only 2 elements, than there is only one event, so no interevent times
    rescaling!(model, params, events) # Transform the time of the events
    X        = diff(events) # Calculate the interevent times)
    sort!(X)
    LIMIT    = 8.0
    STEP     = 0.01
    N        = Float64(length(X))
    integral = 0.0
    @inbounds for x in 0.0:STEP:LIMIT # The integral in approximated as a square function
        integral += exp(-x) *
                    ((sum(exp.(-X .* x)) / N) - # The empirical Laplace transform
                    (1.0 / (1.0 + x)))^2.0 # The Laplace transform of the exponential distribution (λ / (s + λ))
    end
    return integral * STEP
end

#--------------- Laplace distance 2 -------------------------------------------
function statistic(model::ParametricTemporalModel, ::Type{LPExponential2}, params::Parameters, events::Times)
    (length(events) <= 2) && return 1.0 # If 'events' has only 2 elements, than there is only one event, so no interevent times
    X   = diff(rescaling!(model, params, events)[1])
    X ./= mean(X)
    XX  = X .+ X'
    N   = Float64(length(X))
    return (sum((1 .+ ((XX' .+ 2).^2)) ./ ((XX .+ 1).^3)) / N) -
           (2 * sum((X .+ 2) ./ ((X .+ 1).^2))) +
           N
end

#--------------- KS uniform distance ---------------------------------------------
function statistic(model::ParametricTemporalModel, ::Type{KSUniform}, params::Parameters, events::Times)
    (length(events) <= 2) && return 1.0 # If 'events' has only 2 elements, than there is only one event, so no interevent times
    T    = rescaling!(model, params, events)
    N    = Float64(length(events))
    xs   = diff(events)
    vals = ecdf(xs)
    d    = 0.0
    @inbounds for i in eachindex(xs)
        d = max(d, (xs[i] / T) - vals[i] , vals[i+1] - (xs[i] / T))
    end
    return d
end

#--------------- ecdf --------------------------------------------------------
function ecdf(X::AbstractVector{<:Real})
    sort!(xs)
    X = unique!(X)
    N = length(X)
    vals = zeros(length(X) + 1)
    @inbounds for i in eachindex(X)
        vals[i+1] = count(xs .<= X[i]) / N # Calculate the empirical CDF
    end
    return vals
end


