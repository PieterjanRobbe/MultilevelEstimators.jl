## distributions.jl : collection of probability density functions
#
# Representation of some commonly used distributions.
#
# Note: ideally, we want to link with Distributions.jl. However, this is troublesome when
# coupling with QMC points. Minimal requirements are the implementation of rng_native_52 and
# rand, see [RandomNumbers.jl](https://github.com/sunoru/RandomNumbers.jl) and the Julia
# documentation [Random Numbers](https://docs.julialang.org/en/v1/stdlib/Random/index.html).
# As a temporary fix, we will define some commonly used probability distributions, and
# define the required transformations ourselves (with less efficiency).
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## AbstractDistribution ##
"""
AbstractDistribution{T<:AbstractFloat}

Represents an abstract Distribution type that generates values of type `T`.
"""
abstract type AbstractDistribution{T<:AbstractFloat} end

## Uniform ##
"""
```julia
Uniform([a=0., [b=1.]])
``

Returns a uniform distribution on [`a`, `b`].

See also: [`Normal`](@ref), [`TruncatedNormal`](@ref), [`Weibull`](@ref)
"""
struct Uniform{T} <: AbstractDistribution{T}
    a::T
    b::T
end

Uniform(a::T1=0., b::T2=1.) where {T1<:Real, T2<:Real} =
check_finite(a, Uniform, "a") &&
check_finite(b, Uniform, "b") &&
check_ordered(a, b, Uniform, "a", "b") &&
Uniform{promote_type(T1, T2)}(a, b)

show(io::IO, distr::Uniform{T}) where T = print(io, string("Uniform{", T, "}(a=", distr.a, ", b=", distr.b, ")"))

transform(T::Uniform, x::Number) = T.a + (T.b - T.a)*x

## Normal ##
"""
```julia
Normal([μ=0., [σ=1.]])
```

Returns a normal distribution with mean `μ` and standard deviation `σ`.

See also: [`Uniform`](@ref), [`TruncatedNormal`](@ref), [`Weibull`](@ref)
"""
struct Normal{T} <: AbstractDistribution{T}
    μ::T
    σ::T
end

Normal(μ::T1=0., σ::T2=1.) where {T1<:Real, T2<:Real} =
check_finite(μ, Normal, "μ") &&
check_finite(σ, Normal, "σ") &&
check_positive(σ, Normal, "σ") &&
Normal{promote_type(T1, T2)}(μ, σ)

show(io::IO, distr::Normal{T}) where T = print(io, string("Normal{", T, "}(μ=", distr.μ, ", σ=", distr.σ, ")"))

transform(T::Normal, x::Number) = T.μ + T.σ*Φ⁻¹(x)

## TruncatedNormal ##
"""
```julia
TruncatedNormal([μ=0., [σ=1., [a=-2., [b=2.]]]])
```

Returns a truncated normal distribution with mean `μ` and standard deviation `σ`.

See also: [`Uniform`](@ref), [`Normal`](@ref), [`Weibull`](@ref)
"""
struct TruncatedNormal{T} <: AbstractDistribution{T}
    μ::T
    σ::T
    a::T
    b::T
end

TruncatedNormal(μ::T1=0., σ::T2=1., a::T3=-2., b::T4=2.) where {T1<:Real, T2<:Real, T3<:Real, T4<:Real} =
check_finite(μ, TruncatedNormal, "μ") &&
check_finite(σ, TruncatedNormal, "σ") &&
check_finite(a, TruncatedNormal, "a") &&
check_finite(b, TruncatedNormal, "b") &&
check_positive(σ, TruncatedNormal, "σ") &&
check_ordered(a, b, TruncatedNormal, "a", "b") &&
TruncatedNormal{promote_type(T1, T2, T3, T4)}(μ, σ, a, b)

show(io::IO, distr::TruncatedNormal{T}) where T = print(io, string("TruncatedNormal{", T, "}(μ=", distr.μ, ", σ=", distr.σ, ", a=", distr.a, ", b=", distr.b, ")"))

function transform(T::TruncatedNormal, x::Number)
    α = (T.a - T.μ)/T.σ
    β = (T.b - T.μ)/T.σ
    T.μ + T.σ*Φ⁻¹(Φ(α) + (Φ(β) - Φ(α))*x)
end

Φ⁻¹(x::T where {T<:Real}) = √2*erfinv(2*x-1)

Φ(x::T where {T<:Real}) = 1/2*(1+erf(x/√2))

## Weibull ##
"""
```julia
Weibull([k=2., [λ=√2]])
```

Returns a Weibull distribution with shape parameter `k` and scale parameter `λ`.

See also: [`Uniform`](@ref), [`Normal`](@ref), [`TruncatedNormal`](@ref)
"""
struct Weibull{T} <: AbstractDistribution{T}
    k::T
    λ::T
end

Weibull(k::T1=2., λ::T2=√2) where {T1<:Real, T2<:Real} =
check_finite(k, Weibull, "k") &&
check_finite(λ, Weibull, "λ") &&
check_positive(k, Weibull, "k") &&
check_positive(λ, Weibull, "λ") &&
Weibull{promote_type(T1, T2)}(k, λ)

show(io::IO, distr::Weibull{T}) where T = print(io, string("Weibull{", T, "}(k=", distr.k, ", λ=", distr.λ, ")"))

transform(T::Weibull, x::Number) = T.λ * ( -log(1 - x))^(1/T.k)

## checks ##
check_finite(x, type_name, parameter_name) = isfinite(x) || throw(ArgumentError(string("in ", type_name, ", ", parameter_name, " must be finite")))
check_ordered(a, b, type_name, parameter_name_1, parameter_name_2) = a < b || throw(ArgumentError(string("in ", type_name, ", ", parameter_name_1, " must be smaller than ", parameter_name_2)))
check_positive(x, type_name, parameter_name) = x > 0 || throw(ArgumentError(string("in ", type_name, ", ", parameter_name, " must be positive")))
