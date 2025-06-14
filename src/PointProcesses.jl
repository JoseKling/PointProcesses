module PointProcesses

using Random
using Optim
using StaticArrays
using HCubature
using Distributions
using Base.Threads
import HypothesisTests: ExactOneSampleKSTest, pvalue
import StatsBase: sample, mean

include("Processes.jl")
include("Models.jl")
include("TPP.jl")
include("SPP.jl")
include("DomainFunction.jl")
include("ParametricModels.jl")
include("Statistics.jl")
include("FitTests.jl")
include("Tests.jl")

include.(["TemporalModels/" * file for file in readdir(joinpath(@__DIR__, "TemporalModels"))])
include.(["SpatialModels/" * file for file in readdir(joinpath(@__DIR__, "SpatialModels"))])

end
