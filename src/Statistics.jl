# export Raw, Scaling, Thinning, Voronoi

########################## Statistic - data type ###############################
abstract type Statistic end
# struct Scaling  <: Statistic end
# struct Inverse  <: Statistic end
# struct Thinning <: Statistic end
# struct Voronoi  <: Statistic end
# struct Raw      <: Statistic end

# TODO: Pearson residual, Voronoi residual
#---------------- Scaling residual ---------------------------------------
function statistic(model::ParametricModel, ::Type{Scaling}, pp::PointProcess
                   ; n_splits=1000)
    M = 0.0
    ∫λ = ∫CIF(model, pp)
    for i in 1:n_splits
        point = i * measure(pp.A) / n_splits
        n = sum(x .<= point)
        s = abs(n - ∫λ(pp.A[1], point))
        M = s > M ? s : M
    end
    return M
end

#---------------- Inverse residual ---------------------------------------
function statistic(model::ParametricModel, ::Type{Inverse}, pp::PointProcess
                   ; n_splits=100)
    M = 0.0
    λ = CIF(model, pp)
    for i in 1:n_splits
        point = i * measure(pp.A) / n_splits
        s = abs(sum(sqrt.(λ.((@view pp.events[pp.events .<= point])))) - (point - pp.A[1]))
        M = s > M ? s : M
    end
    return M
end

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

function statistic(model::ParametricModel, ::Type{Raw}, pp::PointProcess)
    cif = intensity(model, pp)
    S   = 0.0
    for i in 1:length(pp)
        expected = cif.integral((0.0, pp[i]))
        S        = max(S, expected - (i - 1), i - expected)
    end
    return S
end
