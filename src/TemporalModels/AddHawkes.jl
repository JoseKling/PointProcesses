export AddHawkes

struct AddHawkes <: ParametricModel
    I::Interval
    f
    ∫f
    ∫If::Real

    function AddHawkes(I::Tuple{Real, Real}, f, ∫f=nothing, ∫If=nothing)
        I, ∫f, ∫If = initialize(I, f, ∫f, ∫If)
        new(I, f, ∫f, ∫If)
    end

    AddHawkes(a::Real, b::Real, f, ∫f=nothing, ∫If=nothing) = 
        AddHawkes((a, b), f, ∫f, ∫If)
end

function simulate(model::AddHawkes, params::Parameters{4})
    sim_base = [simulate(params[1], model.I);
                thinning(simulate(params[2], model.I), model.f)]
    return sort!([sim_base;
                 generate_descendants(sim_base,
                                      model.I.b,
                                      params[3],
                                      params[4])])
end

estimate(model::AddHawkes, events::Events) = estimate(model, events, zeros(length(events)), zeros(length(events)))

function estimate(model::AddHawkes, events::Events, lambda_events::Events, f_times::Events) 
    N          = length(events)
    N == 0 && return (0.0, 0.0, 0.0, 0.0)
    T          = model.I.b[1] - model.I.a[1]
    events   .*= (N / T) # Average inter-event time equal to 1. Numerical stability
    step_tol   = 1e-4
    max_iter   = 500
    c1         = 0.2 + (0.6 * rand())
    c2         = 0.1 + (0.8 * rand())
    f_times   .= model.f.(events) # Calculate the f at event times
    μ, γ, ψ, β = (c1 / 2, c1 / 2, (1.0 - c1), log(10.0) * c2) # ψ is the branching factor, α = ψ * β
    iter       = 0
    error      = step_tol + 1.0 # Initial error is larger than tolerance
    while (error < step_tol) && (iter < max_iter) # Stop iteration when optimization is smaller then tolerance
        iter  = iter + 1
        lambda_events[1] = 0.0
        @inbounds for i in 2:N
            lambda_events[i] = exp(β * (events[i] - events[i-1])) * (1 + lambda_events[i-1])
        end
        @. lambda_events = (μ + (γ * f_times)) + (ψ * β * lambda_events)
        D      = 0.0 # Expected number of descendants
        I      = 0.0 # Expected number of immigrants
        I_base = 0.0 # Expected number of immigrants from the base intensity
        div    = 0.0 # For calculating the new parameters for the next iteration
        @inbounds for i in 2:N # Go through each event calculating the probabilities of being an immigrant or a descendant
            D_i  = 0.0 # Probability that t_i is a descendant
            @inbounds for j in 1:i-1
                diffs = events[i] - events[j]
                D_ij  = (ψ * β * exp(β * diffs)) / lambda_events[i] # Probability t_i is a descendant of t_j
                div  += (diffs * D_ij)
                D_i  += D_ij
            end
            D      += D_i
            I      += 1 - D_i
            I_base += (1 - D_i) * (μ / (μ + (γ * f_times[i])))
        end
        error = max(μ - (D / N), γ - (((I - I_base) * T) / (model.∫If * N)), ψ - (D / div), β - (D / div))
        μ, γ, ψ, β = I_base / N, ((I - I_base) * T) / (model.∫If * N), D / N, D / div # Update parameters for the next iteration
    end
    iter >= max_iter && @warn("Maximum number of iterations reached without convergence.") 
    events .*= (T / N) # Unnormalize events to original scale
    return (μ * (N / T),
            γ * (N / T),
            ψ * β * (N / T),
            β * (N / T))
end

function time_transform(model::AddHawkes, params::Parameters{4}, events::Events)
    return time_transform(HH(I), (params.μ, params.α, params.β), events) .+
           time_transform(IP(I, model.f, model.∫f, model.∫If), (0.0, params.γ), events)
end

function CIF(model::AddHawkes, params::Parameters{4}, events::Events)
    sort!(events)
    function cif(t)
        activation = t < events ? 0 : sum(exp.(params[3] .* (@view events[events .< t])))
        return params[1] + (params[2] * model.f(t)) + (params[3] * activation / exp(params[4] * t))
    end
    return cif
end

# function ∫CIF(params::ParametersIH,
#                 rec::Record,
#                 f::Interp)
#     ∫base   = (params.μ * span(rec))
#     ∫f      = (params.γ * ∫(f))
#     ∫self   = (params.α / params.β) * sum(1 .- exp.(-params.β .* (rec.finish .- rec.events)))
#     return ∫base + +∫f + ∫self
# end

function min_max_CIF(model::AddHawkes, params::Parameters{4}, events::Events)
    M = 0
    for i in eachindex(events)
        temp = sum(exp.(-params[4] * (events[i] - events[1:i])))
        M = (temp > M) ? temp : M
    end
    return params[1], params[1] + params[2] + (params[3] * M)
end



# function simulate!(model::AddHawkes, params::Parameters{4}, events::Events)
#     N_base = simulate!(IP(model.I, model.f, model.∫f, model.∫If), (params[1], params[2]), events)
#     N_desc = generate_descendants!(events, N_base, model.I.b, params[3], params[4])
#     return N_desc
# end
