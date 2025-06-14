export KSExponential, LPExponential, KSUniform

struct KSExponential <: Statistic end
struct KSUniform     <: Statistic end
struct LPExponential <: Statistic end

# TODO: Reduce allocations
#---------------- Kolmogorov-Smirnov distance ---------------------------------
function statistic(model::ParametricModel, ::Type{KSExponential}, tpp::TPP)::Float64
    (length(tpp.times) <= 2) && return 1.0 # If 'events' has only 2 elements, than there is only one event, so no interevent times
    X = diff(rescale(model, tpp).times)
    return KSDistance(X, x -> 1 - exp(-x)) # T is the mean of the interevent times
end

#--------------- Laplace distance 2 -------------------------------------------
function statistic(model::ParametricModel, ::Type{LPExponential}, tpp::TPP)::Float64
    (length(tpp.times) <= 2) && return 1.0 # If 'events' has only 2 elements, than there is only one event, so no interevent times
    X   = diff(rescale(model, tpp).times)
    X ./= mean(X)
    XX  = X .+ X'
    N   = Float64(length(X))
    return (sum((1 .+ ((XX' .+ 2).^2)) ./ ((XX .+ 1).^3)) / N) -
           (2 * sum((X .+ 2) ./ ((X .+ 1).^2))) +
           N
end

#--------------- KS uniform distance ---------------------------------------------
function statistic(model::ParametricModel, ::Type{KSUniform}, tpp::TPP)::Float64
    (length(tpp.times) <= 2) && return 1.0 # If 'events' has only 2 elements, than there is only one event, so no interevent times
    rescaled = rescale(model, tpp)
    return KSDistance(diff(rescaled.times), x -> x / rescaled.I[2]) # Calculate the empirical CDF
end

#---------------- Kolmogorov-Smirnov distance calculation -----------------------
function KSDistance(X::Vector{Float64}, f)
    N   = Float64(length(X))
    i   = 1
    val = 0.0
    d   = 0.0
    while true
        lower = f(X[i]) - val
        while i < N && X[i+1] == X[i]
            i += 1
        end
        val = i / N
        upper = val - f(X[i])
        d = max(d, lower, upper)
        i += 1
        i >= N && break
    end
    return d

end
