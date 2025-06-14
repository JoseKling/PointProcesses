#---------------------- Simulation ---------------------------------------
# TODO: Documentation
"""
Simulate a Poisson process. Base algorithm for all other simulations.
"""
function simulate(mu::Real, A::Area{D}) where {D}
    n = rand(Poisson(mu * measure(A)))
    return [Tuple(rand(D) .* (A[2] .- A[1]) .+ A[1]) for _ in 1:n]
end

"""
Implements the thinning algorithm.
[Lewis (1979)](https://doi.org/10.1002/nav.3800260304)
"""
function thinning(locs::Locations, f)
    return locs[rand(length(locs)) .<= f.(locs)] # Simplified for a normalized CIF
end

function thinning(locs::Locations, f, maxf::Real, minf::Real)
    Ui = (rand(length(locs)) .* (maxf - minf)) .+ minf
    return locs[Ui .<= f.(locs)]
end
