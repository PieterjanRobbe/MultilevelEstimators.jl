# TODO: fix printing of \'e
# TODO: evaluation of covariance function is crucial and bottleneck!!!
# TODO: complete evaluation of exponential/squared exponential
# TODO: maybe split Matern, Exponential, Squared Exponential across multiple files?

## covariance_functions.jl : GRF covariance functions

## CovarianceFunction ##
"""
`struct CovarianceFunction{S,T}`

Implements a covariance function of type `T`. The first parameter, `S`, can be `Separable` or `NonSeparable`.
"""
struct CovarianceFunction{S,T} 
    cov::T
end

struct Separable end
struct NonSeparable end

show(io::IO,::Type{Separable}) = print(io,"separable")
show(io::IO,::Type{NonSeparable}) = print(io,"")
show(io::IO, C::CovarianceFunction{S,T}) where {S,T} = print(io, "$(S) $(C.cov)")

## Mat\'ern ##
struct Matern{d,T}
    λ::T
    ν::T
    p::T
end

function Matern(d::N where {N<:Integer},λ::T where {T<:Real},ν::T where {T<:Real},p::T where {T<: Real})
    d > 0 || throw(ArgumentError("dimension must be positive, got $(d)"))
    λ > 0 || throw(ArgumentError("correlation length λ of Mat\'ern covariance cannot be negative or zero!"))
    ν > 0 || throw(ArgumentError("smoothness ν of Mat\'ern covariance cannot be negative or zero!"))
    p >= 1 || throw(ArgumentError("in p-norm, p must be greater than or equal to 1!"))
    Matern{d,promote_type(typeof(λ),typeof(ν),typeof(p))}(promote(λ,ν,p)...) 
end

"""
`struct MaternCovarianceFunction(d, λ, ν, p)`

Create a Mat\'ern covariance function in `d` dimensions with correlation length `λ`, smoothness `ν` and p-norm `p`.

Examples:
```
matern = MaternCovarianceFunction(2, 0.1, 0.5, 2)
julia> CovarianceFunction{NonSeparable,Matern{2,Float64}}(Matern{2,Float64}(0.1, 0.5, 2.0))

```
"""
function MaternCovarianceFunction(args...)
    m = Matern(args...)
    CovarianceFunction{NonSeparable,typeof(m)}(m)
end

"""
`struct SeparableMaternCovarianceFunction(d, λ, ν, p)`

Create a separable Mat\'ern covariance function in `d` dimensions with correlation length `λ`, smoothness `ν` and p-norm `p`.

Examples:
```
matern = SeparableMaternCovarianceFunction(2, 0.1, 0.5, 2)
julia> CovarianceFunction{Separable,Matern{2,Float64}}(Matern{2,Float64}(0.1, 0.5, 2.0))

```
"""
function SeparableMaternCovarianceFunction(args...)
    m = Matern(args...)
    CovarianceFunction{Separable,typeof(m)}(m)
end

function apply(cov::CovarianceFunction{S,Matern} where {S}, x::AbstractMatrix{T}, y::AbstractMatrix{T}) where {T}
    C = zeros(T,size(x,1),size(y,1))
    for i in 1:length(x)
        for j in 1:length(y)
            @inbounds C[i,j] = all(x[i,:].==y[j,:]) ? 1 : 2^(1-cov.ν)/gamma(cov.ν)*(sqrt(2*cov.ν)*norm(x[i]-y[j],cov.p)/cov.λ).^cov.ν.*besselk(cov.ν,sqrt(2*cov.ν)*norm(x[i]-y[j],cov.p)/cov.λ)
        end
    end

    return C
end

show(io::IO,m::Matern{d}) where {d} = print(io, "Mat\'ern covariance function in $(d) dimension"*(d == 1 ? "" : "s")*" with correlation length λ = $(m.λ), smoothness ν = $(m.ν) and $(m.p)-norm")


## Exponential ##
struct Exponential{d,T}
    λ::T
    p::T
end
    
function Exponential(d::N where {N<:Integer},λ::T where {T<:Real},p::T where {T<: Real}) 
    d > 0 || throw(ArgumentError("dimension must be positive, got $(d)"))
    λ > 0 || throw(ArgumentError("correlation length λ of exponential covariance cannot be negative or zero!"))
    p >= 1 || throw(ArgumentError("in p-norm, p must be greater than or equal to 1!"))
    Exponential{d,promote_type(typeof(λ),typeof(p))}(promote(λ,p)...) 
end

"""
`struct ExponentialCovarianceFunction(d, λ, p)`

Create an exponential covariance function in `d` dimensions with correlation length `λ` and p-norm `p`.

Examples:
```
exponential = ExponentialCovarianceFunction(2, 0.1, 2)
julia> CovarianceFunction{NonSeparable,Exponential{2,Float64}}(Exponential{2,Float64}(0.1, 2.0))

```
"""
function ExponentialCovarianceFunction(args...)
    e = Exponential(args...)
    CovarianceFunction{NonSeparable,typeof(e)}(e)
end

"""
`struct SeparableExponentialCovarianceFunction(d, λ, p)`

Create a separable exponential covariance function in `d` dimensions with correlation length `λ` and p-norm `p`.

Examples:
```
exponential = SeparableExponentialCovarianceFunction(2, 0.1, 2)
julia> CovarianceFunction{Separable,Exponential{2,Float64}}(Exponential{2,Float64}(0.1, 2.0))

```
"""
function SeparableExponentialCovarianceFunction(args...)
    e = Exponential(args...)
    CovarianceFunction{Separable,typeof(e)}(e)
end

function apply(cov::CovarianceFunction{S,Exponential} where {S}, x::AbstractMatrix{T}, y::AbstractMatrix{T}) where {T}
    # TODO
end

show(io::IO,e::Exponential{d}) where {d} = print(io, "exponential covariance function in $(d) dimension"*(d == 1 ? "" : "s")*" with correlation length λ = $(e.λ) and $(e.p)-norm")

## SquaredExponential ##
struct SquaredExponential{d,T}
    λ::T
    p::T
end
    
function SquaredExponential(d::N where {N<:Integer},λ::T where {T<:Real},p::T where {T<: Real}) 
    d > 0 || throw(ArgumentError("dimension must be positive, got $(d)"))
    λ > 0 || throw(ArgumentError("correlation length λ of exponential covariance cannot be negative or zero!"))
    p >= 1 || throw(ArgumentError("in p-norm, p must be greater than or equal to 1!"))
    SquaredExponential{d,promote_type(typeof(λ),typeof(p))}(promote(λ,p)...) 
end

"""
`struct SquaredExponentialCovarianceFunction(d, λ, p)`

Create a squared exponential covariance function in `d` dimensions with correlation length `λ` and p-norm `p`.

Examples:
```
sq_exponential = SquaredExponentialCovarianceFunction(2, 0.1, 2)
julia> CovarianceFunction{NonSeparable,SquaredExponential{2,Float64}}(SquaredExponential{2,Float64}(0.1, 2.0))

```
"""
function SquaredExponentialCovarianceFunction(args...)
    e = SquaredExponential(args...)
    CovarianceFunction{NonSeparable,typeof(e)}(e)
end

"""
`struct SeparableSquaredExponentialCovarianceFunction(d, λ, p)`

Create a separable squared exponential covariance function in `d` dimensions with correlation length `λ` and p-norm `p`.

Examples:
```
sq_exponential = SeparableSquaredExponentialCovarianceFunction(2, 0.1, 2)
julia> CovarianceFunction{Separable,SquaredExponential{2,Float64}}(SquaredExponential{2,Float64}(0.1, 2.0))

```
"""
function SeparableSquaredExponentialCovarianceFunction(args...)
    e = SquaredExponential(args...)
    CovarianceFunction{Separable,typeof(e)}(e)
end

function apply(cov::CovarianceFunction{S,SquaredExponential} where {S}, x::AbstractMatrix{T}, y::AbstractMatrix{T}) where {T}
    # TODO
end

show(io::IO,e::SquaredExponential{d}) where {d} = print(io, "squared exponential covariance function in $(d) dimension"*(d == 1 ? "" : "s")*" with correlation length λ = $(e.λ) and $(e.p)-norm")
