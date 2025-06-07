#---------------------- Simulation ---------------------------------------
# TODO: Documentation
"""
Simulate a Poisson process. Base algorithm for all other simulations.
"""
function simulate(mu::Real, A::Area{d}) where {d}
    n = rand(Poisson(mu * measure(A)))
    return [Tuple(rand(d) .* (A[2] .- A[1]) .+ A[1]) for _ in 1:n]
end

"""
Implements the thinning algorithm.
[Lewis (1979)](https://doi.org/10.1002/nav.3800260304)
"""
function thinning(events::Locations, f)
    return events[rand(length(events)) .<= f.(events)] # Simplified for a normalized CIF
end

function thinning(events::Locations, f, maxf::Real, minf::Real)
    Ui = (rand(length(events)) .* (maxf - minf)) .+ minf
    return events[Ui .<= f.(events)] # Simplified for a normalized CIF
end
