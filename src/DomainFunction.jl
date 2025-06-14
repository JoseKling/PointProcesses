export DomainFunction

struct DomainFunction{T<:Union{Float64, Tuple{Vararg{Float64}}}, F, IF}
    f::F
    domain::NTuple{2, T}
    integral::IF

    function DomainFunction(f::F, domain::NTuple{2, T}; integral=nothing) where {F, T}
        new_domain = (Float64.(domain[1]), Float64.(domain[2]))  # Ensure domain is in Float64
        new_T = typeof(new_domain[1])
        !(f(new_domain[1]) isa Float64) && throw(TypeError(:DomainFunction, "f, return type of function", Float64, typeof(f(new_domain[1]))))
        if integral === nothing
            default_integral = _default_integral(f, new_T)
            return new{new_T, F, typeof(default_integral)}(f, new_domain, default_integral)
        else
            !(integral(new_domain) isa Float64) && throw(TypeError(:DomainFunction, "integral, return type of integral", Float64, typeof(f(new_domain[1]))))
            return new{new_T, F, typeof(integral)}(f, new_domain, integral)
        end
    end
end
(df::DomainFunction)(x) = df.f(x)

# Default integral for 1D (temporal)
function _default_integral(f, ::Type{Float64})
    (I; rtol=1e-6) -> hquadrature(f, I[1], I[2]; rtol=rtol)[1]
end

# Default integral for ND (spatial)
function _default_integral(f, ::Type{NTuple{D, Float64}}) where D
    return (A; rtol=1e-6) -> hcubature(f, A[1], A[2]; rtol=rtol)[1]
end