## distribution.jl : collection of probability density functions
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
# This file is part of MultilevelEstimators.jl - A Julia toolbox for
# Multilevel Monte Carlo Methods (c) Pieterjan Robbe, 2019

## AbstractDistribution ##
abstract type AbstractDistribution end

## Uniform ##
"""
    Uniform([a = 0, [b = 1]])

Return a uniform distribution on [`a`, `b`].

# Examples
```jldoctest; setup = :(using MultilevelEstimators)
julia> Uniform()
Uniform{Int64}(a = 0, b = 1)
```
See also: [`Normal`](@ref), [`TruncatedNormal`](@ref), [`Weibull`](@ref)
"""
struct Uniform{T} <: AbstractDistribution
    a::T
    b::T
end

Uniform(a::T1=0, b::T2=1) where {T1<:Real, T2<:Real} =
check_finite(Uniform, a, "a") &&
check_finite(Uniform, b, "b") &&
check_ordered(Uniform, a, b, "a", "b") &&
Uniform{promote_type(T1, T2)}(a, b)

show(io::IO, distr::Uniform{T}) where T = print(io, string("Uniform{", T, "}(a = ", distr.a, ", b = ", distr.b, ")"))

"""
    transform(D::AbstractDistribution, x::Real)

Apply a transformation to the uniformly distributed number `x` such that it is distributed according to `D`.

The transform is implemented by applying the inverse CDF of the distribution `D`.

!!! note

    This is currently limited only to distributions with an analytic expression for the inverse CDF.

# Examples
```jldoctest; setup = :(using MultilevelEstimators; import Random; Random.seed!(1))
julia> x = rand(10)
10-element Array{Float64,1}:
 0.23603334566204692
 0.34651701419196046
 0.3127069683360675
 0.00790928339056074
 0.4886128300795012
 0.21096820215853596
 0.951916339835734
 0.9999046588986136
 0.25166218303197185
 0.9866663668987996

julia> broadcast(i -> transform(Normal(), i), x)
10-element Array{Float64,1}:
 -0.7191204773101185 
 -0.39474101736576256
 -0.4881918784059056 
 -2.413074939355834  
 -0.02854727903738184
 -0.8030663180302484 
  1.6637253880181684 
  3.7310515426407296 
 -0.6692682660790877 
  2.2163540181904553 
```
"""
transform(T::Uniform, x::Number) = T.a + (T.b - T.a)*x

## Normal ##
"""
    Normal([μ = 0, [σ = 1]])

Return a normal distribution with mean `μ` and standard deviation `σ`.

# Examples
```jldoctest; setup = :(using MultilevelEstimators)
julia> Normal()
Normal{Int64}(μ = 0, σ = 1)
```
See also: [`Uniform`](@ref), [`TruncatedNormal`](@ref), [`Weibull`](@ref)
"""
struct Normal{T} <: AbstractDistribution
    μ::T
    σ::T
end

Normal(μ::T1=0, σ::T2=1) where {T1<:Real, T2<:Real} =
check_finite(Normal, μ, "μ") &&
check_finite(Normal, σ, "σ") &&
check_larger_than(Normal, σ, "σ", 0) &&
Normal{promote_type(T1, T2)}(μ, σ)

show(io::IO, distr::Normal{T}) where T = print(io, string("Normal{", T, "}(μ = ", distr.μ, ", σ = ", distr.σ, ")"))

transform(T::Normal, x::Number) = T.μ + T.σ*Φ⁻¹(x)

## TruncatedNormal ##
"""
    TruncatedNormal([μ = 0, [σ = 1, [a = -2, [b = 2]]]])

Return a truncated normal distribution with mean `μ` and standard deviation `σ`.

# Examples
```jldoctest; setup = :(using MultilevelEstimators)
julia> TruncatedNormal()
TruncatedNormal{Int64}(μ = 0, σ = 1, a = -2, b = 2)
```
See also: [`Uniform`](@ref), [`Normal`](@ref), [`Weibull`](@ref)
"""
struct TruncatedNormal{T} <: AbstractDistribution
    μ::T
    σ::T
    a::T
    b::T
end

TruncatedNormal(μ::T1=0, σ::T2=1, a::T3=-2, b::T4=2) where {T1<:Real, T2<:Real, T3<:Real, T4<:Real} =
check_finite(TruncatedNormal, μ, "μ") &&
check_finite(TruncatedNormal, σ, "σ") &&
check_finite(TruncatedNormal, a, "a") &&
check_finite(TruncatedNormal, b, "b") &&
check_larger_than(TruncatedNormal, σ, "σ", 0) &&
check_ordered(TruncatedNormal, a, b, "a", "b") &&
TruncatedNormal{promote_type(T1, T2, T3, T4)}(μ, σ, a, b)

show(io::IO, distr::TruncatedNormal{T}) where T = print(io, string("TruncatedNormal{", T, "}(μ = ", distr.μ, ", σ = ", distr.σ, ", a = ", distr.a, ", b = ", distr.b, ")"))

function transform(T::TruncatedNormal, x::Number)
    α = (T.a - T.μ)/T.σ
    β = (T.b - T.μ)/T.σ
    T.μ + T.σ*Φ⁻¹(Φ(α) + (Φ(β) - Φ(α))*x)
end

Φ⁻¹(x::T where {T<:Real}) = √2*erfinv(2*x-1)

Φ(x::T where {T<:Real}) = 1/2*(1+erf(x/√2))

## Weibull ##
"""
    Weibull([k = 2, [λ = √2]])

Return a Weibull distribution with shape parameter `k` and scale parameter `λ`.

# Examples
```jldoctest; setup = :(using MultilevelEstimators)
julia> Weibull()
Weibull{Float64}(k = 2.0, λ = 1.4142135623730951)
```
See also: [`Uniform`](@ref), [`Normal`](@ref), [`TruncatedNormal`](@ref)
"""
struct Weibull{T} <: AbstractDistribution
    k::T
    λ::T
end

Weibull(k::T1=2, λ::T2=√2) where {T1<:Real, T2<:Real} =
check_finite(Weibull, k, "k") &&
check_finite(Weibull, λ, "λ") &&
check_larger_than(Weibull, k, "k", 0) &&
check_larger_than(Weibull, λ, "λ", 0) &&
Weibull{promote_type(T1, T2)}(k, λ)

show(io::IO, distr::Weibull{T}) where T = print(io, string("Weibull{", T, "}(k = ", distr.k, ", λ = ", distr.λ, ")"))

transform(T::Weibull, x::Number) = T.λ * ( -log(1 - x))^(1/T.k)
