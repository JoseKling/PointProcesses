abstract type Model           end
abstract type ParametricModel end

abstract type ParametricTemporalModel          <: ParametricModel end
abstract type ParametricSpatialModel{d}        <: ParametricModel end
abstract type ParametricMarkModel{T}           <: ParametricModel end
abstract type ParametricSpatioTemporalModel{d} <: ParametricModel end

#------------ Model initialization ----------------------
convert_interval(a::Real, b::Real; d::Int=0) = convert_interval((a, b))
convert_interval(a::Tuple{Vararg{<:Real}}, b::Tuple{Vararg{<:Real}}; d::Int=0) = convert_interval((a, b))
function convert_interval(; d::Int=0)
    d == 0 && return (0.0, 1.0)
    return (Tuple(zeros(d)), Tuple(ones(d)))
end
function convert_interval(I::Tuple{Real, Real})
    if I[1] > I[2]
        return Float64.((I[2], I[1]))
    else
        return Float64.(I)
    end
end
function convert_interval(A::NTuple{2, Tuple{Vararg{<:Real}}})
    mins = map(min, A[1], A[2])
    maxs = map(max, A[1], A[2])
    return (Float64.(mins), Float64.(maxs))
end

initialize(a::Real, b::Real, args...; d::Int=0) =
    initialize(convert_interval(a, b), args...; d=d)
initialize(a::NTuple{2, Real}, b::NTuple{2, Real}, args...; d::Int=0) =
    initialize(convert_interval(a, b), args...; d=d)
initialize(f, ∫f=nothing; d::Int=0) = initialize(convert_interval(d=d), f, ∫f)
function initialize(I::Union{Interval, Area}, f, ∫f=nothing; d::Int=0)
    I   = convert_interval(I)
    @assert f(I[1]) isa Float64 "Function must return a Float64."
    xs = random_points(I, 100000)
    if !isapprox(maximum(f.(xs)), 1.0, atol=1e-2) || !isapprox(minimum(f.(xs)), 0.0, atol=1e-2)
        @warn "Methods assume f i normalized to min(f) = 0 and max(f) = 1. Are you sure this holds?"
    end
    if isnothing(∫f)
        @warn "Not providing the integral of f may cause performance and precision problems"
        ∫f = integral(f)
        ∫If = ∫f(I, 1000000)
    else
        intf = integral(f)
        max_diff = maximum([∫f((I[1], xs[i])) - intf((I[1], xs[i])) for i in 1:50])
        max_diff > 1e-1 && @warn "Is the integral of f correct?"
        ∫If = ∫f(I)
    end
    return I, f, ∫f, ∫If
end

function integral(f)
    return (A, n=1000) -> mapreduce(f, +, random_points(A, n)) * measure(A) / n
end
random_points(I::Interval, n::Int) = (rand(n) * measure(I)) .+ I[1]
random_points(A::Area{d}, n::Int) where d = Tuple.(eachcol(A[1] .+ (rand(d, n) .* (A[2] .- A[1]))))

