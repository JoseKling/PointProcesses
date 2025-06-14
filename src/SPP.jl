export SPP

const Locations{D, T} = Vector{NTuple{D, T}} where T<:Real
const Area{D, T} = NTuple{2, NTuple{D, T}} where T<:Real
measure(A::Area{D, T}) where {D, T<:Real} = prod(abs.(A[2] .- A[1]))

struct SPP{D} <: SpatialProcess
    locations::Locations{D, Float64}
    region::Area{D, Float64}

    function SPP(locations::Locations{D, Float64}, region::Area{D, Float64}) where D
        new_region = new_area(region)
        new_locs = new_locations(locations, new_region)
        return new{D}(new_locs, new_region)
    end
end
Base.getindex(spp::SPP, i::Int) = spp.locations[i]

function new_area(region::Area{D, Float64}) where D
    all(region[1] .== region[2]) && error("Area must have non-zero size.")
    mins = map(min, region[1], region[2])
    maxs = map(max, region[1], region[2])
    return (mins, maxs)
end

function new_locations(locations::Locations{D, Float64}, region::Area{D, Float64}) where D
    if any([any(x .< region[1]) || any(x .> region[2]) for x in locations])
        @warn "Some events are outside the area [$(region[1]), $(region[2])]. These events will be ignored."
        new_locs = filter(x -> all(x .>= region[1]) && all(x .<= region[2]), locations)
    else
        new_locs = locations
    end
    return sort(new_locs, by=x->x)
end
