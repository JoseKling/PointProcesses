abstract type PointProcess end

abstract type TemporalProcess <: PointProcess end
abstract type SpatialProcess  <: PointProcess end
# TODO: Add MarkedProcess and SpatioTemporalProcess types
# abstract type MarkedProcess   <: PointProcess end
# abstract type SpatioTemporalProcess <: PointProcess end

Base.length(pp::PointProcess) = length(getfield(pp, 1))
Base.getindex(pp::PointProcess, I::AbstractVector{<:Integer}) = [pp[i] for i in I]
Base.show(io::IO, pp::PointProcess) = print(io, "Point process over $(pp.region) with $(length(pp)) events.")