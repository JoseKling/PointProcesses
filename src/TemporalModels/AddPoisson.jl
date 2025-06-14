function simulate(model::AddPoisson, region::Interval)
    times = simulate(model.μ, region)
    append!(times, thinning(simulate(model.γ, region), model.f))
    sort!(times)
    return TPP(times, region)
end

function estimate(model::AddPoisson, tpp::TPP)
    isempty(tpp.times) && return AddPoisson((0.0, 0.0))
    f_times   = model.f.(tpp.times)
    T         = measure(tpp.region)
    initial_x = [length(f_times) / T, length(f_times) / T]
    lower     = [0.0, 0.0]
    upper     = [Inf, Inf]
    ∫If       = integral(model.f, tpp.region)
    results   = optimize(x -> objective_ap(x, f_times, ∫If, T),
                         (stor, x) -> gradient_ap!(stor, x, f_times, ∫If, T),
                         (stor, x) -> hessian_ap!(stor, x, f_times),
                         lower,
                         upper,
                         initial_x,
                         IPNewton())
    return AddPoisson(Optim.minimizer(results)[1], Optim.minimizer(results)[2], model.f)
end

function objective_ap(x, f_times, ∫f, M)
    return (x[1] .* M) + (x[2] .* ∫f) - # integral of the CIF 
           sum(log.(x[1] .+ (x[2] .* f_times))) # Sum of log of f at event times
end

function gradient_ap!(storage, x, f_times, ∫f, M)
    storage[1] = M - sum(1.0 ./ (x[1] .+ (x[2] .* f_times)))
    storage[2] = ∫f - sum(f_times ./ (x[1] .+ (x[2] .* f_times)))
end

function hessian_ap!(storage, x, f_times)
    storage[1, 1] = sum(1.0 ./ ((x[1] .+ (x[2] .* f_times))) .^ 2.0)
    storage[1, 2] = sum(f_times ./ ((x[1] .+ (x[2] .* f_times))) .^ 2.0)
    storage[2, 1] = storage[1, 2]
    storage[2, 2] = sum((f_times .^ 2) ./ ((x[1] .+ (x[2] .* f_times))) .^ 2.0)
end

function rescale(model::AddPoisson, tpp::TPP)
    times = model.μ .* tpp.times
    @. times .+= [model.γ * model.f.integral.(tpp.region[1], t) for t in times]
    return TPP(times, (0.0, (measure(tpp.region) * model.μ) + (model.γ * model.f.integral.(tpp.region))))
end

function intensity(model::AddPoisson, tpp::TPP)
    cif(t) = model.μ + (model.γ * model.f(t)), tpp.region
    ∫cif(region) = (model.μ * measure(region)) + (model.γ * model.f.integral(region))
    return DomainFunction(cif, tpp.region, integral=∫cif)
end

function simulate!(model::AddPoisson, tpp::TPP, sim::Times)
    n_const = simulate!(params[1], tpp.region, sim)
    n_f     = simulate!(params[2], tpp.region, (@view sim[n_const+1: end]))
    n_f     = thinning!((@view sim[n_const+1: n_const+n_f]), model.f)
    sort!((@view sim[1: n_const+n_f]))
    return n_const + n_f
end

function rescale!(model::AddPoisson, tpp::TPP, rescaled::Times)
    @. rescaled[1:length(tpp.times)] = model.μ * tpp.times
    for i in eachindex(tpp.times)
        rescaled[i] += model.γ * model.f.∫f((0.0, tpp.times[i]))
    end
end