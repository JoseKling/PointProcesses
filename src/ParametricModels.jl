################ Interface for Parametric Models ################
export simulate, estimate, rescale, intensity

function simulate(::ParametricModel, region::Union{Interval, Area})
    error("Simulation not implemented for $(typeof(model))")
end

function estimate(::ParametricModel, pp::PointProcess)
    error("Estimation not implemented for $(typeof(model))")
end

function intensity(::ParametricModel, pp::PointProcess)
    error("Intensity not implemented for $(typeof(model))")
end

function rescale(model::ParametricModel, pp::PointProcess)
    error("rescale not implemented for $(typeof(model))")
end

############### Parametric Models ###############
export HPoisson, HHawkes, AddPoisson, AddHawkes, MultPoisson

struct HPoisson <: ParametricModel
    μ::Float64

    HPoisson() = new(NaN)
    HPoisson(μ) = new(μ)
end

struct HHawkes <: ParametricModel
    μ::Float64
    α::Float64
    β::Float64

    HHawkes() = new(NaN, NaN, NaN)
    HHawkes(μ, α, β) = new(μ, α, β)
end

struct AddPoisson <: ParametricModel
    μ::Float64
    γ::Float64
    f::DomainFunction

    AddPoisson(f::DomainFunction) = new(NaN, NaN, f)
    AddPoisson(μ, α, f::DomainFunction) = new(μ, α, f)
end

struct AddHawkes <: ParametricModel where T
    μ::Float64
    γ::Float64
    α::Float64
    β::Float64
    f::DomainFunction

    AddHawkes(f::DomainFunction) = new(NaN, NaN, NaN, NaN, f)
    AddHawkes(μ, γ, α, β, f::DomainFunction) = new(μ, γ, α, β, f)
end

struct MultPoisson <: ParametricModel
    μ::Float64
    γ::Float64
    f::DomainFunction

    MultPoisson(f::DomainFunction) = new(NaN, NaN, f)
    MultPoisson(μ, α, f::DomainFunction) = new(μ, α, f)
end