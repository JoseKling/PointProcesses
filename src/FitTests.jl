export fit_test, no_bootstrap
export reynaud2014

#---------------------- Goodness-of-fit test --------------------
function fit_test(model::Model, S::Type{<:Statistic}, pp::PointProcess; n_sims::Integer=1000)
    model_est = estimate(model, pp)
    s         = statistic(model_est, S, pp)
    sim_stats = Vector{Float64}(undef, n_sims)
    Threads.@threads for i in 1:n_sims
        sim       = simulate(model_est, pp.A)
        model_sim = estimate(model, sim)
        sim_stats[i] = statistic(model_sim, S, sim)
    end
    ns = count(>=(s), sim_stats)
    return (ns + 1) / (n_sims + 1)
end

function no_bootstrap(model::Model, S::Type{<:Statistic}, pp::PointProcess; n_sims::Integer=1000)
    model_est = estimate(model, pp)
    s         = statistic(model_est, S, pp)
    sim_stats = Vector{Float64}(undef, n_sims)
    Threads.@threads for i in 1:n_sims
        sim          = simulate(model_est, pp.A)
        sim_stats[i] = statistic(model_est, S, sim)
    end
    ns = count(>=(s), sim_stats)
    return (ns + 1) / (n_sims + 1) # The p-value of the test
end

function reynaud2014(model::ParametricModel, tpp::TPP)
    params_est = estimate(model, tpp)
    rescaled   = rescale(params_est, tpp)
    n_S        = round(Int, (length(rescaled.times) - 1)^(2/3))
    S          = sample(diff(rescaled.times), n_S, replace=false)
    return pvalue(ExactOneSampleKSTest(S, Exponential()))
end
