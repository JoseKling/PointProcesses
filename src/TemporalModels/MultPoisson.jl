# TODO: Bug in simulation
function simulate(model::MultPoisson, region::Interval)
    mu_ip    = model.γ < 0.0 ? model.μ * exp(model.γ) : model.μ
    gamma_ip = model.γ < 0.0 ? model.μ - mu_ip : model.μ * (exp(model.γ) - 1.0)
    sim = simulate(mu_ip, region)
    append!(sim, thinning(simulate(gamma_ip, region), CIF(model, (mu_ip, mu_ip + gamma_ip))))
    sort!(sim)
    return TPP(sim, region)
end

function estimate(model::MultPoisson, tpp::TPP)
    isempty(tpp.times) && return (0.0, 0.0)
    f_times = model.f.(tpp.times)
    divs = 1000
    r = rand()
    f_mesh = map(model.f, mesh(tpp.region, divs))
    m = measure(tpp.region) / divs
    initial_x = [r, 1 - r]
    lower     = [0.0, -Inf]
    upper     = [Inf, Inf]
    results   = optimize(x -> objective_mp(x, f_times, f_mesh, m),
                         (stor, x) -> gradient_mp!(stor, x, f_times, f_mesh, m),
                         (stor, x) -> hessian_mp!(stor, x, f_times, f_mesh, m),
                         lower,
                         upper,
                         initial_x,
                         IPNewton())
    return MultPoisson(Optim.minimizer(results)[1], Optim.minimizer(results)[2], model.f)
end

function objective_mp(x, f_times, f_mesh, m)
    return (x[1] * (m * sum(exp.(x[2] .* f_mesh)))) - # integral of the CIF 
           sum(log.(x[1] .* exp.(x[2] .* f_times))) # Sum of log of f at event times
end

function gradient_mp!(storage, x, f_times, f_mesh, m)
    storage[1] = (m * sum(exp.(x[2] .* f_mesh))) - (length(f_times) / x[1])
    storage[2] = (x[1] * (m * sum(f_mesh .* exp.(x[2] .* f_mesh)))) - sum(f_times)
end

function hessian_mp!(storage, x, f_times, f_mesh, m)
    storage[1, 1] = length(f_times) / x[1]^2
    storage[1, 2] = m * sum(f_mesh .* exp.(x[2] .* f_mesh))
    storage[2, 1] = storage[1, 2]
    storage[2, 2] = x[1] * (m * sum((f_mesh.^2) .* exp.(x[2] .* f_mesh)))
end

intensity(model::MultPoisson, tpp::TPP) = DomainFunction(x -> model.μ * exp(model.γ * model.f(x)), tpp.region)
