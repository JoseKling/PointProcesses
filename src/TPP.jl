export TPP

const Times{T} = Vector{T} where T<:Real
const Interval{T} = NTuple{2, T} where T<:Real
measure(region::Interval) = abs(region[2] - region[1])

struct TPP <: TemporalProcess
    times::Times{Float64}
    region::Interval{Float64}

    function TPP(times::Times, region::Interval)
        new_region = new_interval(region)
        new_ts = new_times(times, new_region)
        return new(new_ts, new_region)
    end
end
Base.getindex(tpp::TPP, i::Int) = tpp.times[i]

function new_interval(I::Interval)
    I[1] == I[2] && error("Interval must have non-zero length.")
    return I[1] > I[2] ? (I[2], I[1]) : I
end

function new_times(times::Times, I::Interval)
    new_times = times
    if any([x < I[1] || x > I[2] for x in times])
        @warn "Some events are outside the interval [$(I[1]), $(I[2])]. These events will be ignored."
        new_times = filter(x -> x >= I[1] && x <= I[2], times)
    end
    return sort(new_times)
end

function interevents(tpp::TPP)
    diffs = similar(tpp.times)
    diffs[1] = tpp.times[1] - tpp.I[1]  # First event relative to the start of the interval
    @inbounds for i in 2:length(tpp.times)
        diffs[i] = tpp.times[i] - tpp.times[i-1]
    end
    return diffs
end

function rescale(model::TemporalModel, tpp::TPP)
    cif = intensity(model, tpp)
    return TPP([cif.integral(tpp.I[1], t) for t in tpp.times],
               (0.0, cif.integral(tpp.I)))
end