function simulate(model::AddPoisson, region::Area)
    return [simulate(model.μ, region); thinning(simulate(model.γ, region), model.f)]
end

function estimate(model::AddPoisson, spp::SPP)
    isempty(spp.locations) && return AddPoisson(0.0, 0.0, model.f)
    f_times   = model.f.(spp.locations)
    initial_x = [length(f_times) / measure(spp.region), length(f_times) / measure(spp.region)]
    lower     = [0.0, 0.0]
    upper     = [Inf, Inf]
    results   = optimize(x -> objective_sap(x, f_times, model.∫Af, measure(model.region)),
                         (stor, x) -> gradient_sap!(stor, x, f_times, model.∫If, measure(model.region)),
                         (stor, x) -> hessian_sap!(stor, x, f_times),
                         lower,
                         upper,
                         initial_x,
                         IPNewton())
    return AddPoisson(Optim.minimizer(results)[1], Optim.minimizer(results)[2], model.f)
end

function objective_sap(x, f_times, ∫f, M)
    return (x[1] .* M) + (x[2] .* ∫f) - # integral of the CIF 
           sum(log.(x[1] .+ (x[2] .* f_times))) # Sum of log of f at event times
end

function gradient_sap!(storage, x, f_times, ∫f, M)
    storage[1] = M - sum(1.0 ./ (x[1] .+ (x[2] .* f_times)))
    storage[2] = ∫f - sum(f_times ./ (x[1] .+ (x[2] .* f_times)))
end

function hessian_sap!(storage, x, f_times)
    storage[1, 1] = sum(1.0 ./ ((x[1] .+ (x[2] .* f_times))) .^ 2.0)
    storage[1, 2] = sum(f_times ./ ((x[1] .+ (x[2] .* f_times))) .^ 2.0)
    storage[2, 1] = storage[1, 2]
    storage[2, 2] = sum((f_times .^ 2) ./ ((x[1] .+ (x[2] .* f_times))) .^ 2.0)
end

function intensity(model::AddPoisson, spp::SPP)
    cif(x) = model.μ + (model.γ * model.f(x))
    ∫cif(region) = (model.μ * measure(region)) + (model.γ * model.∫f(region))
    return DomainFunction(cif, spp.region, integral=∫cif)
end
