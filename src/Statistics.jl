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


function statistic(model::ParametricModel, ::Type{Raw}, params::Parameters, events::Events)
    cif = ∫CIF(model, params)
    S   = 0.0
    for i in eachindex(events)
        expected = cif((0.0, events[i]))
        S        = max(S, expected - (i - 1), i - expected)
    end
    return S
end
