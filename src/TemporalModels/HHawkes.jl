function simulate(model::HHawkes, region::Interval)
    sim = simulate(model.μ, region)
    sim_desc = generate_descendants(sim, region[2], model.α, model.β)
    append!(sim, sim_desc)
    sort!(sim)
    return TPP(sim, region)
end

function estimate(::HHawkes, tpp::TPP; step_tol::Float64 = 1e-3, max_iter::Int = 1000)
    N         = length(tpp.times)
    N        == 0 && return HHawkes((0.0, 0.0, 0.0))
    T         = measure(tpp.region)
    c1        = 0.2 + (0.6 * rand())
    c2        = 0.1 + (0.8 * rand())
    norm_ts   = tpp.times .* (N / T) # Average inter-event time equal to 1. Numerical stability
    μ, ψ, β   = c1, (1.0 - c1), log(10.0) * c2 # ψ is the branching factor, α = ψ * β
    iter      = 0
    lambda_events = zeros(N) # λ(t_i) = μ + ψ * β * exp(-β * (t_i - t_{i-1}))
    error = step_tol + 1.0 # Initial error is larger than tolerance
    while (error >= step_tol) #&& (iter < max_iter) # Stop iteration when optimization is smaller then tolerance 
        iter += 1
        lambda_events[1] = 0.0
        @inbounds for i in 2:N
            lambda_events[i] = exp(-β * (norm_ts[i] - norm_ts[i-1])) * (1.0 + lambda_events[i-1])
        end
        lambda_events .*= (ψ * β)
        lambda_events .+= μ
        D   = 0.0 # Expected number of descendants
        div = 0.0 # Needed to calculate the new parameters for the next iteration
        @inbounds for i ∈ 2:N
            @inbounds for j ∈ 1:i-1
                diffs = norm_ts[i] - norm_ts[j]
                D_ij  = (ψ * β * exp(-β * diffs)) / lambda_events[i] # Probability that t_i is a descendant of t_j
                D    += D_ij
                div  += (diffs * D_ij)
            end
        end
        error = max(abs(μ - (1 - (D / N))), abs(ψ - (D / N)), abs(β - (D / div)))
        μ, ψ, β = 1 - (D / N),  D / N, D / div # Update parameters for the next iteration
    end
    iter >= max_iter && @warn("Maximum number of iterations reached without convergence.")
    return HHawkes(μ * (N / T), ψ * β * (N / T), β * (N / T)) # Unnormalize parameters end
end

function rescale(model::HHawkes, tpp::TPP)
    N = length(tpp.times)
    A = zeros(N + 1) # Array A in Ozaki (1979)
    @inbounds for i in 2:N
        A[i] = exp(-model.β * (tpp.times[i] - tpp.times[i-1])) * (1 + A[i-1])
    end
    A[end] = exp(-model.β * (tpp.region[2] - tpp.times[end])) * (1 + A[end-1])
    times  = model.μ .* (tpp.times .- tpp.region[1]) # Integral of base rate
    T_base = model.μ * measure(tpp.region)
    @inbounds for i in eachindex(times)
        times[i] += (model.α / model.β) * ((i - 1) - A[i])
    end
    return TPP(times, (0.0, T_base + ((model.α / model.β) * (N - A[end]))))
end

# TODO: Implement the integral
function intensity(model::HHawkes, tpp::TPP)
    function cif(t)
        activation = sum(exp.(model.β .* (@view tpp.times[1:searchsortedfirst(tpp.times, t) - 1])))
        return model.μ + (model.α * activation / exp(model.β * t))
    end
    return DomainFunction(cif, tpp.region)
end

function simulate!(model::HHawkes, region::Interval, sim::Vector{Float64})
    N_base = simulate!(model.μ, region, sim)
    N_desc = generate_descendants!(sim, N_base, region[2], model.α, model.β)
    sort!((@view sim[1:N_base + N_desc]))
    return N_base + N_desc
end
