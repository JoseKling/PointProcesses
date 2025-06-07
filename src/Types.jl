export Events, Locations, Interval, Area, Parameters
export KSExponential, LPExponential, LPExponential2, KSUniform
export Raw, Scaling, Thinning, Voronoi

########################## Basic data types #################################
const Time  = Float64
const Times = AbstractVector{Float64}
const Interval = Tuple{Float64, Float64}

const Location{d} = NTuple{d, Float64}
const Locations{d} = Vector{Location{d}}
const Area{d} = Tuple{Location{d}, Location{d}}

const Events = Union{Times, Locations}

const Parameters{d} = NTuple{d, Float64}

########################## Interval - data type ###############################

measure(I::Interval) = I[2] - I[1]
measure(A::Area) = prod(A[2] .- A[1])


########################## Statistic - data type ###############################
abstract type Statistic end
struct Scaling       <: Statistic end
struct Inverse       <: Statistic end
struct Pearson       <: Statistic end
struct Thinning      <: Statistic end
struct Voronoi       <: Statistic end
struct Raw           <: Statistic end
struct KSExponential <: Statistic end
struct KSUniform     <: Statistic end
struct LPExponential <: Statistic end
struct LPExponential2 <: Statistic end

