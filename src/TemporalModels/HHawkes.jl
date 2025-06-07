export HHawkes

struct HHawkes <: ParametricTemporalModel
    I::Interval

    HHawkes(arg...) = new(convert_interval(arg...))
end

function simulate(model::HHawkes, params::Parameters{3})
    sim_μ = simulate(params[1], model.I)
    return sort!([sim_μ;
                 generate_descendants(sim_μ,
                                      model.I[2],
                                      params[2],
                                      params[3])])
end

function estimate(model::HHawkes, events::Events; step_tol::Float64 = 1e-3, max_iter::Int = 1000)
    N         = length(events)
    N        == 0 && return (0.0, 0.0, 0.0)
    T         = model.I[2] - model.I[1]
    c1        = 0.2 + (0.6 * rand())
    c2        = 0.1 + (0.8 * rand())
    events  .*= (N / T) # Average inter-event time equal to 1. Numerical stability
    μ, ψ, β   = c1, (1.0 - c1), log(10.0) * c2 # ψ is the branching factor, α = ψ * β
    iter      = 0
    lambda_events = zeros(N) # λ(t_i) = μ + ψ * β * exp(-β * (t_i - t_{i-1}))
    exp_events    = zeros(N)
    error = step_tol + 1.0 # Initial error is larger than tolerance
    while (error >= step_tol) #&& (iter < max_iter) # Stop iteration when optimization is smaller then tolerance 
        iter += 1
        lambda_events[1] = 0.0
        @inbounds for i in 2:N
            lambda_events[i] = exp(-β * (events[i] - events[i-1])) * (1.0 + lambda_events[i-1])
        end
        lambda_events .*= (ψ * β)
        lambda_events .+= μ
        D   = 0.0 # Expected number of descendants
        div = 0.0 # Needed to calculate the new parameters for the next iteration
        @inbounds for i ∈ 2:N
            @inbounds for j ∈ 1:i-1
                diffs = events[i] - events[j]
                D_ij  = (ψ * β * exp(-β * diffs)) / lambda_events[i] # Probability that t_i is a descendant of t_j
                D    += D_ij
                div  += (diffs * D_ij)
            end
        end
        error = max(abs(μ - (1 - (D / N))), abs(ψ - (D / N)), abs(β - (D / div)))
        μ, ψ, β = 1 - (D / N),  D / N, D / div # Update parameters for the next iteration
    end
    iter >= max_iter && @warn("Maximum number of iterations reached without convergence.")
    events .*= (T / N) # Unnormalize events to original scale
    return (μ * (N / T),
            ψ * β * (N / T),
            β * (N / T)) # Unnormalize parameters
end

function rescaling!(model::HHawkes, params::Parameters{3}, events::Events)
    N         = length(events)
    A         = zeros(N + 1) # Array A in Ozaki (1979)
    @inbounds for i in 2:N
        A[i] = exp(-params[3] * (events[i] - events[i-1])) * (1 + A[i-1])
    end
    A[end]  = exp(-params[3] * (model.I[2] - events[end])) * (1 + A[end-1])
    events, T_base .= rescaling(params[1], events, model.I) # Integral of base rate
    @inbounds for i in eachindex(events)
        events[i] += (params[2] / params[3]) * ((i - 1) - A[i])
    end
    return T_base + ((params[2] / params[3]) * (N - A[end]))
end

function CIF(model::HHawkes, params::Parameters{3}, events::Events)
    sort!(events)
    function cif(t)
        @assert model.I[1] <= t <= model.I[2] "Time $t is outside the interval $(model.I[1]), $(model.I[2])"
        activation = t < events[1] ? 0 : sum(exp.(params[3] .* (@view events[events .< t])))
        return params[1] + (params[2] * activation / exp(params[3] * t))
    end
    return cif
end

