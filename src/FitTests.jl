export fit_test, no_bootstrap
export reynaud2014

#---------------------- Goodness-of-fit test --------------------
"""
Performs the goodness-of-fit test for a given model and distance function. The 
number of simulations used in the bootstrap can be set with the `n_sims` keyword
argument (default 1000).

The `model` may be provided either as a string or as an instance of the model type.
The `dist` may be provided either as a string or as an instance of the distance type.

For inhomogeneous processes, the `proxy` field of `rec` must not be `nothing`.

Returns a named tuple with fields:
- `p`         -> returned p-value of the test
- `sim_dists` -> simulated distances used in calculating the p-value
- `dist`      -> distance between the observed and the estimated process
- `params`    -> estimated parameters of the model
See [`Record`](@ref) [`distance`](@ref) [`Model`](@ref) [`simulate`](@ref) [`estimate`](@ref)

```@example
rec = Record("Resources/data/record.csv")
fit_test("hp", "ks", rec) # Perform the goodness-of-fit test for a Poisson process using the KS-distance
fit_test("hh", "lp", rec; n_sims=1000) # Perform the test for a Hawkes process using the Laplace distance.
```
"""
function fit_test(model::ParametricModel, S::Type{<:Statistic}, events::Events; n_sims::Integer=1000)
    params = estimate(model, events) # Estimate parameters via MLE
    s      = statistic(model, S, params, events)
    s_sim  = zeros(n_sims)
    for i in eachindex(s_sim) # For each of the simulated processes
        sim        = simulate(model, params)
        sim_params = estimate(model, sim) # Estimate parameters via MLE
        s_sim[i]   = statistic(model, S, sim_params, sim)
    end
    p_value = (sum(s_sim .>= s) + 1) / (n_sims + 1)
    return (p         = p_value,
            params    = params,
            stats     = s,
            sim_stats = s_sim)
end

function no_bootstrap(model::ParametricModel, S::Type{<:Statistic}, events::Events; n_sims::Integer=1000)
    params = estimate(model, events) # Estimate parameters via MLE
    s      = statistic(model, S, params, events)
    s_sim  = zeros(n_sims)
    for i in eachindex(s_sim) # For each of the simulated processes
        sim      = simulate(model, params)
        s_sim[i] = statistic(model, S, params, sim)
    end
    p_value = (sum(s_sim .>= s) + 1) / (n_sims + 1)
    return (p         = p_value,
            params    = params,
            stats     = s,
            sim_stats = s_sim)
end

function reynaud2014(model::ParametricModel, events::Events)
    params    = estimate(model, events) # Estimate parameters via MLE
    transf, _ = rescaling(model, params, events)
    n_S       = round(Int, (length(transf) - 1)^(2/3))
    S         = sample(diff(transf), n_S, replace=false)
    return pvalue(ExactOneSampleKSTest(S, Exponential()))
end
