#----------- rescaling ------------------
function rescaling!(model::ParametricTemporalModel, params::Parameters, events::Times)
    @warn "rescaling not implemented. Using standard implementation."
    ∫cif = ∫CIF(model, params, events)
    events .= [∫cif((0.0, event)) for event in events]
    return ∫cif(model.I)
end

"""
    rescaling(model::ParametricModel, params::Parameters, events::Times)

Performs the time re-scaling \$t \\mapsto \\int_0^t \\lambda(s) ds\$, where \$\\lambda(s)\$ is the
conditional intensity function. The resulting process is a homogeneous Poisson process.
The function returns a vector containing the transformed events and the measure of the re-scaled interval.
"""
function rescaling(model::ParametricTemporalModel, params::Parameters{n}, events::Times) where n
    transformed = zeros(length(events))
    T = rescaling!(transformed)
    return transformed, T
end

"""
    rescaling(μ::Real, events::Times, I::Interval)

Time re-scaling for a homogeneous Poisson process with rate `μ` over the interval `I`.
Used in multiple models.
"""
function rescaling(μ::Real, events::Times, I::Interval)
    rescaled = zeros(length(events))
    T = rescaling!(rescaled)
    return rescaled, T
end

"""
    rescaling!(μ::Real, events::Times, I::Interval)

Time re-scaling for a homogeneous Poisson process with rate `μ` over the interval `I`.
Used in multiple models. In-place version that modifies `events`.
"""
function rescaling!(μ::Real, events::Times, I::Interval)
    @. events = (events - I[1]) * μ
    return μ * measure(I)
end

# ----------- Simulation ----------------------------
"""
    simulate(mu::Float64, T::Float64)

Simulate a Poisson process over the interval \$[0, T]\$ and intensity mu.
"""
function simulate(mu::Float64, T::Float64)
    n = rand(Poisson(mu * T))
    return sort!(rand(n) .* T)
end

simulate(mu::Float64, I::Interval) = I[1] .+ simulate(mu, measure(I))

"""
    thinning(events::Times, f)

Thinning algorithm for a normalized function \$f\$, in which
\$\\min(f) = 0\$ and \$\\max(f) = 1\$.
[Lewis (1979)](https://doi.org/10.1002/nav.3800260304)

"""
function thinning(events::Times, f)
    return events[rand(length(events)) .<= f.(events)] # Simplified for a normalized CIF
end

"""
    thinning(events::Times, f)

Thinning algorithm for a function \$f\$ with \$\\min(f) =\$minf and
\$\\max(f) =\$minf.
[Lewis (1979)](https://doi.org/10.1002/nav.3800260304)
"""
function thinning(events::Times, f, maxf::Real, minf::Real)
    Ui = (rand(length(events)) .* (maxf - minf)) .+ minf
    return events[Ui .<= f.(events)] # Simplified for a normalized CIF
end

function simulate_inverse(∫If::Float64, inv∫f)
    sim = simulate(1.0, ∫If)
    return inv∫f.(sim)
end

function simulate_single_inverse(inv∫f)
    t_star = randexp()
    return inv∫f(t_star)
end

"""
Generates the descendants of a Hawkes process from the events generated from the
baseline intensity (the `parents`).

The idea is to recurse over the generations. The offspring of the current `parents`
will become the next parents until some generation does not produce any offspring.

For each event in `parents`, new events are generated as an inhomogeneous Poisson
process with the activation function.
"""
function generate_descendants(parents::Times,
                              T::Float64,
                              α::Float64,
                              β::Float64)
    descendants  = eltype(parents)[]
    next_gen     = parents
    a_over_b     = α / β
    b_over_a     = β / α
    b_inv        = 1.0 / β
    while ~isempty(next_gen)
        last_gen = copy(next_gen)
        next_gen = eltype(parents)[]
        for parent in last_gen
            invf(t) = parent - (b_inv * log(1.0 - (b_over_a * t)))
            sim_transf = simulate_inverse(a_over_b * (1.0 - exp(β * (parent - T))), invf)
            append!(next_gen, sim_transf)
        end
        append!(descendants, next_gen)
    end
    return descendants
end












# """
#     simulate!(model::ParametricModel, params::Parameters, sim::Events)

# Simulates events from the model given the parameters and stores them in `sim`.
# The type of `params` should match the model's requirements. The function returns
# the number of simulated events.
# """
# function simulate!(model::ParametricModel, params::Parameters{n}, sim::Events) where n
#     error("Simulation not implemented for $(typeof(model))")
# end

# function generate_descendants!(events::Events,
#                                N::Integer,
#                                T::Float64,
#                                α::Float64,
#                                β::Float64)
#     pos_gen = 1
#     N_gen   = N
#     b_inv        = -1.0 / β
#     b_over_a     = β / α
#     a_over_b     = α / β
#     while N_gen != 0
#         pos_desc = pos_gen + N_gen
#         for parent in (@view events[pos_gen:pos_gen + N_gen - 1])
#             N_desc = simulate!(1.0,
#                                a_over_b * (1.0 - exp(β * (parent - T))),
#                                (@view events[pos_desc:end])) # Simulate a homogeneous Poisson process
#             @inbounds for i in pos_desc:pos_desc + N_desc - 1
#                 events[i] = parent - (b_inv * log(1.0 - (b_over_a * events[i])))
#             end
#             pos_desc = pos_desc + N_desc
#         end
#         pos_gen = pos_gen + N_gen
#         N_gen   = pos_desc - pos_gen
#     end
#     sort!(@view events[1:pos_gen - 1])
#     return pos_gen - 1
# end

# function thinning!(events::Events, f)
#     n = 1
#     for event in events
#         if rand() <= f(event)
#             events[n] = event
#             n += 1
#         end
#     end
#     return n - 1
# end


# function simulate!(mu::Float64, T::Float64, sim::Events)
#     n = 1
#     mu_inv = 1 / mu
#     t = -log(rand()) * mu_inv
#     while t <= T
#         sim[n] = t
#         n += 1
#         t += -log(rand()) * mu_inv
#     end
#     return n - 1
# end

# function simulate!(mu::Real, I::Interval, sim::Events)
#     n = simulate!(mu, (I[2] - I[1]), sim)
#     (@view sim[1:n]) .+= I[1]
#     return n
# end

