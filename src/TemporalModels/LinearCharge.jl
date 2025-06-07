export LinearCharge

struct LinearCharge <: ParametricTemporalModel
    I::Interval

    LinearCharge(arg...) = new(convert_interval(arg...))
end

function simulate(model::LinearCharge, params::Parameters{3})
    n      = 1
    a      = params[3] / 2
    events = Float64[]
    λt_i   = Float64[]
    t_star = randexp()
    t      = (-params[1] + sqrt((params[1]^2) + (2 * params[3] * t_star))) / params[3]
    λ      = params[1] + (params[3] * t)
    while t < model.I.b
        push!(events, t)
        push!(λt_i, λ)
        t_star += randexp()
        b       = params[1] - (params[2] * sum(λt_i))
        c       = (params[2] * sum(λt_i .* events)) - t_star
        t       = (-b + sqrt((b^2) - (4 * a * c))) / (2 * a) # Inverse transform of charge function
        λ       = (λt_i[end] * (1 - params[2])) + (params[3] * (t - events[end]))
        n      += 1
    end
    return events
end

function estimate(model::LinearCharge, events::Times)
    isempty(events) && return (0.0, 0.0)
    initial_x = [length(events) / measure(model.I), rand(), rand() + 1]
    vec       = map(x -> x[1], events)
    lower     = [0.0, 0.0, 0.0]
    upper     = [Inf, 1.0, Inf]
    results   = optimize(x -> objective_charge(x, vec, model.I.a[1], model.I.b[1]),
                         (stor, x) -> gradient_charge!(stor, x, vec, model.I.a[1], model.I.b[1]),
                         (stor, x) -> hessian_charge!(stor, x, vec, model.I.a[1], model.I.b[1]),
                         lower,
                         upper,
                         initial_x,
                         IPNewton())
    return (Optim.minimizer(results)[1],
            Optim.minimizer(results)[2],
            Optim.minimizer(results)[3])
end

function objective_charge(x, ev, a, b)
    λ = zeros(length(ev))
    λ[1] = x[1] .+ (x[3] .* ev[1])
    for i in 2:length(ev)
        λ[i] = (λ[i-1] * (1 - x[2])) + (x[3] * (ev[i] - ev[i-1]))
    end
    return (x[1] * (b - a)) + (x[3] * ((b - a)^2) / 2) - (x[2] * sum(λ .* (b .- ev))) - # integral of the CIF 
    sum(log.(x[1] .+ (x[3] * (ev.^2) / 2) .- (x[2] .* [sum(λ[1:j]) for j in 0:(length(ev)-1)]))) # Sum of log of f at event times
end

function gradient_charge!(storage, x, ev, a, b)
    λ    = zeros(length(ev))
    λ[1] = x[1] .+ (x[3] .* ev[1])
    for i in 2:length(ev)
        λ[i] = (λ[i-1] * (1 - x[2])) + (x[3] * (ev[i] - ev[i-1]))
    end
    sums_λ     = [sum(λ[1:j]) for j in 0:(length(λ)-1)]
    divs       = 1 ./ (x[1]  .+ (x[3] .* ev) .- (x[2] .* sums_λ))
    storage[1] = (b - a) - sum(1.0 .* divs)
    storage[2] = -sum(λ .* (b .- ev)) + sum(sums_λ .* divs)
    storage[3] = (x[3] * ((b - a)^2) / 2) - sum(ev .* divs)
end

function hessian_charge!(storage, x, ev, a, b)
    λ = zeros(length(ev))
    λ[1] = x[1] .+ (x[3] .* ev[1])
    for i in 2:length(ev)
        λ[i] = (λ[i-1] * (1 - x[2])) + (x[3] * (ev[i] - ev[i-1]))
    end
    sums_λ        = [sum(λ[1:j]) for j in 0:(length(ev)-1)]
    divs          = 1 ./ (x[1] .- (x[2] .* sums_λ) .+ (x[3] .* ev)).^2
    storage[1, 1] = sum(1.0 ./ divs)
    storage[2, 2] = -sum((sums_λ.^2) .* divs)
    storage[3, 3] = sum((ev.^2) ./ divs)
    storage[1, 2] = sum(sums_λ .* divs)
    storage[2, 1] = storage[1, 2]
    storage[1, 3] = sum(ev .* divs)
    storage[3, 1] = storage[1, 3]
    storage[2, 3] = sum(sums_λ .* ev .* divs)
    storage[3, 2] = storage[2, 3]
end

# function rescaling!(model::LinearCharge, params::Parameters{3}, events::Times)
# end

function CIF(model::LinearCharge, params::Parameters{3}, events::Times)
    λt_i = zeros(length(events))
    λt_i[1] = params[1] .+ (params[3] .* events[1][1])
    for i in 2:length(events)
        λt_i[i] = (λt_i[i-1] * (1 - params[2])) + (params[3] * (events[i][1] - events[i-1][1]))
    end
    function f(x)
        idx = searchsortedfirst(map(e -> e[1], events), x[1])
        if idx == 1
            return params[1] + (params[3] * x[1])
        else
            return (λt_i[idx-1] * (1 - params[2])) + (params[3] * (x[1] - events[idx-1][1]))
        end
    end
    return f
end

function ∫CIF(model::LinearCharge, params::Parameters{3}, events::Times)
    a   = model.I.a[1]
    b   = model.I.b[1]
    vec = sort(map(x -> x[1], events))
    return x -> (params[1] * measure(model.I)) +
                (params[3] * (b^2 - a^2)) -
                (params[2] * (sum(b .- vec[vec .< b]) - sum(a .- vec[vec .< a])))
end
