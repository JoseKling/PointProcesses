export MultPoisson

struct MultPoisson <: ParametricModel
    I::Interval
    f
    ∫f
    ∫If::Real

    function MultPoisson(I::Tuple{Real, Real}, f, ∫f=nothing, ∫If=nothing)
        I, ∫f, ∫If = initialize(I, f, ∫f, ∫If)
        new(I, f, ∫f, ∫If)
    end

    MultPoisson(a::Real, b::Real, f, ∫f=nothing, ∫If=nothing) = 
        MultPoisson((a, b), f, ∫f, ∫If)
end

function simulate(model::MultPoisson, params::Parameters{2})
    mu_ip    = params[2] < 0.0 ? params[1] * exp(params[2]) : params[1]
    gamma_ip = params[2] < 0.0 ? params[1] - mu_ip : params[1] * (exp(params[2]) - 1.0)
    return sort!([simulate(mu_ip, model.I);
        thinning(simulate(gamma_ip, model.I), CIF(model, params), mu_ip, mu_ip + gamma_ip)])
end

function estimate(model::MultPoisson, events::Events)
    isempty(events) && return (0.0, 0.0)
    f_times = model.f.(events)
    divs = 1000
    r = rand()
    f_mesh = map(model.f, mesh(model.I, divs))
    m = measure(model.I) / divs
    initial_x = [r, 1 - r]
    lower     = [0.0, -Inf]
    upper     = [Inf, Inf]
    results   = optimize(x -> objective_mip(x, f_times, f_mesh, m),
                         (stor, x) -> gradient_mip!(stor, x, f_times, f_mesh, m),
                         (stor, x) -> hessian_mip!(stor, x, f_times, f_mesh, m),
                         lower,
                         upper,
                         initial_x,
                         IPNewton())
    return (Optim.minimizer(results)[1], Optim.minimizer(results)[2])
end

function objective_mip(x, f_times, f_mesh, m)
    return (x[1] * (m * sum(exp.(x[2] .* f_mesh)))) - # integral of the CIF 
           sum(log.(x[1] .* exp.(x[2] .* f_times))) # Sum of log of f at event times
end

function gradient_mip!(storage, x, f_times, f_mesh, m)
    storage[1] = (m * sum(exp.(x[2] .* f_mesh))) - (length(f_times) / x[1])
    storage[2] = (x[1] * (m * sum(f_mesh .* exp.(x[2] .* f_mesh)))) - sum(f_times)
end

function hessian_mip!(storage, x, f_times, f_mesh, m)
    storage[1, 1] = length(f_times) / x[1]^2
    storage[1, 2] = m * sum(f_mesh .* exp.(x[2] .* f_mesh))
    storage[2, 1] = storage[1, 2]
    storage[2, 2] = x[1] * (m * sum((f_mesh.^2) .* exp.(x[2] .* f_mesh)))
end

# function time_transform(model::MultPoisson, params::Parameters{2}, events::Events)
#     hist_ext = [I.a; hist; I.b]
#     xs = LinRange(I.a[1], I.b[1], 1000)
#     ∫cif = integral(Proxy(xs, params.μ .* exp.(params.γ .* f(xs)), normalize=false))
#     transf = ∫cif(hist_ext[2:end])
#     return transf
# end

function CIF(model::MultPoisson, params::Parameters{2}, _...)
    return x -> params[1] * exp(params[2] * model.f(x))
end

function ∫CIF(model::MultPoisson, params::Parameters{2})
    return integral(CIF(model, params), d)
end