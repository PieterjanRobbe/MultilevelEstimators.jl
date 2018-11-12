## distributions.jl : collection of probability density functions
#
# Representation of some commonly used distributions.
#
# Note: ideally, we want to link with Distributions.jl. However, this is troublesome when
# coupling with QMC points. Minimal requirements are the implementation of rng_native_52 and
# rand, see [RandomNumbers.jl](https://github.com/sunoru/RandomNumbers.jl) and the Julia
# documentation [Random Numbers](https://docs.julialang.org/en/v1/stdlib/Random/index.html).
# As a temporary fix, we will define some commonly used probability distributions, and
# define the required transformations ourselves.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## Distribution ##
abstract type Distribution{T<:Number} end

## Uniform ##
struct Uniform{T} <: Distribution{T}
    a::T
    b::T
end

Uniform(a::T1=0, b::T2=1) where {T1, T2} = Uniform{promote_type(T1, T2)}(a, b)

transform(T::Uniform, x::Number) = T.a + (T.b - T.a)*x

## Normal ##
struct Normal{T} <: Distribution{T}
    μ::T
    σ::T
end

Normal(μ::T1=0, σ::T2=1) where {T1, T2} = Normal{promote_type(T1, T2)}(μ, σ)

transform(T::Normal, x::Number) = T.μ + T.σ*Φ⁻¹(x)

## Truncated Normal ##
struct TruncatedNormal{T} <: Distribution{T}
    μ::T
    σ::T
    a::T
    b::T
end

TruncatedNormal(μ::T1=0, σ::T2=1, a::T3=-2, b::T4=2) where {T1, T2, T3, T4} = TruncatedNormal{promote_type(T1, T2, T3, T4)}(μ, σ, a, b)

function transform(T::TruncatedNormal, x::Number)
    α = (T.a - T.μ)/T.σ
    β = (T.b - T.μ)/T.σ
    T.μ + T.σ*Φ⁻¹(Φ(α) + (Φ(β) - Φ(α))*x)
end

Φ⁻¹(x::T where {T<:Real}) = √2*erfinv(2*x-1)

Φ(x::T where {T<:Real}) = 1/2*(1+erf(x/√2))

## Weibull ##
struct Weibull{T} <: Distribution{T}
    k::T
    λ::T
end

Weibull(k::T1=2, λ::T2=√2) where {T1, T2} = Weibull{promote_type(T1, T2)}(k, λ)

transform(T::Weibull, x::Number) = T.λ * ( -log(1 - x))^(1/T.k)
