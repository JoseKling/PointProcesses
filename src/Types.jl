export Event, Events, Parameters, Interval
export KSExponential, LPExponential, KSUniform, Scaling, Thinning, Voronoi
export LPExponential2, Raw

########################## Events - data type #################################
const Event  = Float64
const Events = AbstractVector{Float64}

########################## Parameters - data type #############################
const Parameters{d} = NTuple{d, Float64}

########################## Mesh - data type #################################
struct Mesh
    v1::Event
    v2::Event
    N::Int
end

function (mesh::Mesh)(i::Int) 
    @assert (i >= 1 && i <= mesh.N) "Index out of bounds: $i. Must be between 1 and $(mesh.N)."
    return mesh.v1 + ((2 * i) - 1) * (mesh.v2 - mesh.v1) / (2 * mesh.N)
end

measure(mesh::Mesh) = (mesh.v2 - mesh.v1) / mesh.N
########################## Interval - data type ###############################
Interval = Tuple{Float64, Float64}

measure(I::Interval) = I[2] - I[1]

mesh(I::Interval, N::Int) = mesh(I[1], I[2], N)
function mesh(v1::Event, v2::Event, N::Int)
    return Mesh(v1, v2, N)
end

########################## Statistic - data type ###############################
abstract type Statistic end
struct Scaling       <: Statistic end
struct Inverse       <: Statistic end
struct Pearson       <: Statistic end
struct Thinning      <: Statistic end
struct Voronoi       <: Statistic end
struct KSExponential <: Statistic end
struct KSUniform     <: Statistic end
struct LPExponential <: Statistic end
struct LPExponential2 <: Statistic end
struct Raw <: Statistic end
