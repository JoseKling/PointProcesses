abstract type Model           end
abstract type ParametricModel <: Model end

abstract type TemporalModel          <: Model end
abstract type SpatialModel{d}        <: Model end
abstract type MarkModel{T}           <: Model end
abstract type SpatioTemporalModel{d} <: Model end

abstract type ParametricTemporalModel          <: ParametricModel end
abstract type ParametricSpatialModel{d}        <: ParametricModel end
abstract type ParametricMarkModel{T}           <: ParametricModel end
abstract type ParametricSpatioTemporalModel{d} <: ParametricModel end

const Parameters{d} = NTuple{d, Float64}