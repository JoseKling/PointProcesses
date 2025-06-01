export statistic

#---------------- Scaling residual ---------------------------------------
function statistic(model::ParametricModel, ::Type{Scaling}, params::Parameters, events::Events
                   ; n_splits=1000)
    M = 0.0
    ∫λ = ∫CIF(model, params)
    for i in 1:n_splits
        point = i * measure(model.I) / n_splits
        n = sum(x .<= point)
        s = abs(n - ∫λ(model.I.a, point))
        M = s > M ? s : M
    end
    return M
end

#---------------- Inverse residual ---------------------------------------
function statistic(model::ParametricModel, ::Type{Inverse}, params::Parameters, events::Events
                   ; n_splits=100)
    M = 0.0
    λ = CIF(model, params)
    for i in 1:n_splits
        point = i * measure(I) / n_splits
        s = abs(sum(sqrt.(λ.((@view events[events .<= point])))) - (point - model.I.a))
        M = s > M ? s : M
    end
    return M
end

# TODO: Integral of the square root of the CIF
# #---------------- Pearson residual ---------------------------------------
# function statistic(model::Model,
#                    ::Type{Pearson},
#                    params::Parameters,
#                    events::Events
#                    ; n_splits=100)
#     M = 0.0
#     a = model.I.a
#     b = model.I.b
#     λ = CIF(model, params)
#     for i in 1:n_splits
#         point = i * measure(I) / n_splits
#         s = abs(sum(λ.((@view events[events .<= point]))) - (point - model.I.a))
#         M = s > M ? s : M
#     end
#     return M
# end

#---------------- Thinning residual ---------------------------------------
function statistic(model::ParametricModel, ::Type{Thinning}, params::Parameters, events::Events
                   ; n_thinns=100)
    M = 0.0
    minλ, maxλ = min_max_CIF(model, params, events)
    λ = CIF(model, params, events)
    for i in 1:n_thinns
        thinn = (@view events[rand(length(events)) .<= (minλ ./ λ.(events))])
        N     = length(thinn)
        lower = maximum((thinn ./ measure(model.I)) .- ((0:N-1) / N))
        upper = maximum(((1:N) / N) .- (thinn ./ measure(model.I)))
        s     = max(lower, upper)
        M     = s > M ? s : M
    end
    return M
end

#---------------- Voronoi residual ---------------------------------------
function statistic(model::ParametricModel, ::Type{Voronoi}, params::Parameters, events::Events)
    N = length(events)
    M = 0.0
    ∫λ = ∫CIF(model, params, events)
    voronoi = [model.I.a; (events[1:end-1] .+ events[2:end]) ./ 2; model.I.b]
    integrals = [∫λ((voronoi[i-1], voronoi[i])) for i in 2:length(voronoi)]
    residuals = (1 .- integrals) ./ sqrt.(integrals)
    return maximum(residuals)
end


#---------------- Kolmogorov-Smirnov distance ---------------------------------
function statistic(model::ParametricModel, ::Type{KSExponential}, params::Parameters, events::Events)
    (length(events) <= 2) && return 1.0 # If 'events' has only 2 elements, than there is only one event, so no interevent times
    xs        = diff(time_transform(model, params, events)[1])
    vals = ecdf(xs)
    d        = 0.0
    for i in eachindex(xs)
        exp_cdf = 1.0 - exp.(-xs[i]) # The empirical CDF of the exponential distribution
        d = max(d, exp_cdf - vals[i], vals[i+1] - exp_cdf)
    end
    return d
end

#--------------- Laplace distance ---------------------------------------------
function statistic(model::ParametricModel, ::Type{LPExponential}, params::Parameters, events::Events)
    (length(events) <= 2) && return 1.0 # If 'events' has only 2 elements, than there is only one event, so no interevent times
    time_transform!(model, params, events) # Transform the time of the events
    X        = diff(events) # Calculate the interevent times)
    sort!(X)
    LIMIT    = 8.0
    STEP     = 0.01
    N        = Float64(length(X))
    integral = 0.0
    for x in 0.0:STEP:LIMIT # The integral in approximated as a square function
        integral += exp(-x) *
                    ((sum(exp.(-X .* x)) / N) - # The empirical Laplace transform
                    (1.0 / (1.0 + x)))^2.0 # The Laplace transform of the exponential distribution (λ / (s + λ))
    end
    return integral * STEP
end

#--------------- Laplace distance 2 -------------------------------------------
function statistic(model::ParametricModel, ::Type{LPExponential2}, params::Parameters, events::Events)
    (length(events) <= 2) && return 1.0 # If 'events' has only 2 elements, than there is only one event, so no interevent times
    X   = diff(time_transform(model, params, events)[1])
    X ./= mean(X)
    XX  = X .+ X'
    N   = Float64(length(X))
    return (sum((1 .+ ((XX' .+ 2).^2)) ./ ((XX .+ 1).^3)) / N) -
           (2 * sum((X .+ 2) ./ ((X .+ 1).^2))) +
           N
end

#--------------- KS uniform distance ---------------------------------------------
function statistic(model::ParametricModel, ::Type{KSUniform}, params::Parameters, events::Events)
    (length(events) <= 2) && return 1.0 # If 'events' has only 2 elements, than there is only one event, so no interevent times
    transf, T  = time_transform(model, params, events)
    N          = Float64(length(transf))
    dist_lower = maximum((transf ./ T) .- ((0:N-1) ./ N)) # Largest distance to lower part of step function
    dist_upper = maximum(((1:N) ./ N) .- (transf ./ T)) # Largest distance to upper part of step function
    return max(dist_lower, dist_upper)
end

#--------------- ecdf --------------------------------------------------------
function ecdf(X::AbstractVector{<:Real})
    sort!(xs)
    X = unique!(X)
    N = length(X)
    vals = zeros(length(X) + 1)
    for i in eachindex(X)
        vals[i+1] = count(xs .<= X[i]) / N # Calculate the empirical CDF
    end
    return vals
end

function statistic(model::ParametricModel, ::Type{Raw}, params::Parameters, events::Events)
    cif = ∫CIF(model, params)
    S   = 0.0
    for i in eachindex(events)
        expected = cif((0.0, events[i]))
        S        = max(S, expected - (i - 1), i - expected)
    end
    return S
end