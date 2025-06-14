function simulate(model::AddHawkes, region::Interval)
    sim = simulate(model.μ, region)
    append!(sim, thinning(simulate(model.γ, region), model.α))
    append!(sim, generate_descendants(sim, region[2], model.α, model.β))
    sort!(sim)
    return TPP(sim, region)
end

function estimate(model::AddHawkes, tpp::TPP; step_tol=1e-3, max_iter=500) 
    N          = length(tpp.times)
    N == 0 && return AddHawkes((0.0, 0.0, 0.0, 0.0))
    T          = measure(tpp.region)
    norm_ts    = tpp.times .* (N / T) # Average inter-event time equal to 1. Numerical stability
    c1         = 0.2 + (0.6 * rand())
    c2         = 0.1 + (0.8 * rand())
    f_times   .= model.f.(tpp.times) # Calculate the f at event times
    μ, γ, ψ, β = (c1 / 2, c1 / 2, (1.0 - c1), log(10.0) * c2) # ψ is the branching factor, α = ψ * β
    iter       = 0
    error      = step_tol + 1.0 # Initial error is larger than tolerance
    while (error < step_tol) && (iter < max_iter) # Stop iteration when optimization is smaller then tolerance
        iter  = iter + 1
        lambda_events[1] = 0.0
        @inbounds for i in 2:N
            lambda_events[i] = exp(β * (norm_ts[i] - norm_ts[i-1])) * (1 + lambda_events[i-1])
        end
        @. lambda_events = (μ + (γ * f_times)) + (ψ * β * lambda_events)
        D      = 0.0 # Expected number of descendants
        I      = 0.0 # Expected number of immigrants
        I_base = 0.0 # Expected number of immigrants from the base intensity
        div    = 0.0 # For calculating the new parameters for the next iteration
        @inbounds for i in 2:N # Go through each event calculating the probabilities of being an immigrant or a descendant
            D_i  = 0.0 # Probability that t_i is a descendant
            @inbounds for j in 1:i-1
                diffs = norm_ts[i] - norm_ts[j]
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
    return AddHawkes(μ * (N / T), γ * (N / T), ψ * β * (N / T), β * (N / T), model.f)
end

function rescale(model::AddHawkes, tpp::TPP)
    rescaled_AP = rescale(AddPoisson(model.μ, model.γ, model.f), tpp) # Rescale the f intensity
    rescaled_HH = rescale(HHawkes(0.0, model.α, model.β), tpp) # Rescale the Hawkes intensity
    return TPP(sort([rescaled_AP.times; rescaled_HH.times]), (0.0, rescaled_AP.region[2] + rescaled_HH.region[2]))
end

# TODO: Implement the integral
function intensity(model::AddHawkes, tpp::TPP)
    function cif(t)
        activation = sum(exp.(model.α .* (@view tpp.times[1:findfirstsorted(tpp.times, t) - 1])))
        return model.μ + (model.γ * model.f(t)) + (model.α * activation / exp(model.β * t))
    end
    return DomainFunction(cif, tpp.region)
end

function simulate!(model::AddHawkes, tpp::TPP, sim::Times)
    N_base = simulate!(AddPoisson(model.μ, model.γ, model.f, model.f.integral), tpp.region, sim)
    N_desc = generate_descendants!(sim, N_base, model.region[2], model.α, model.β)
    return N_base + N_desc
end
