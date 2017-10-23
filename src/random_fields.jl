## random_fields.jl : gaussian random fields covariance functions

## abstract covariance functions ##
abstract type CovarianceFunction{d} end

abstract type AbstractMaternCovarianceFunction{d} <: CovarianceFunction{d} end

abstract type AbstractExponentialCovarianceFunction{d} <: CovarianceFunction{d} end

## Matern covariance functions ##
"""
MaternCovarianceFunction

Representation of a Matern covariance function.
"""
struct MaternCovarianceFunction{d,T} <: AbstractMaternCovarianceFunction{d} where T<:Real
    λ::T
    σ::T
    ν::T
    p::T
end

"""
MaternCovarianceFunction(d, λ, σ, ν, p)

Create a Matern covariance function in `d` dimensions with correlation length `λ`, marginal standard deviation `σ`, smoothness `ν` and p-norm `p`.
"""
function MaternCovarianceFunction(d::N where N<:Integer,λ::T where T<:Real,σ::T where T<:Real,ν::T where T<:Real,p::T where T<: Real)  
    d > 0 || throw(ArgumentError("dimension must be positive, got $(d)"))
    λ > 0 || throw(ArgumentError("correlation length λ of the Gaussian random field cannot be negative or zero!"))
    σ > 0 || throw(ArgumentError("marginal standard deviation σ of the Gaussian random field cannot be negative or zero!"))
    ν > 0 || throw(ArgumentError("smoothness ν of Gaussian random field cannot be negative or zero!"))
    p >= 1 || throw(ArgumentError("in p-norm, p must be greater than one!"))

    MaternCovarianceFunction{d,promote_type(typeof(λ),typeof(σ),typeof(ν),typeof(p))}(promote(λ,σ,ν,p)...)
end

"""
SeparableMaternCovarianceFunction

Representation of a separable Matern covariance function.
"""
struct SeparableMaternCovarianceFunction{d,T} <: AbstractMaternCovarianceFunction{d} where T<:Real
    λ::T
    σ::T
    ν::T
    p::T
end

"""
SeparableMaternCovarianceFunction(d, λ, σ, ν, p)

Create a separable Matern covariance function in `d` dimensions with correlation length `λ`, marginal standard deviation `σ`, smoothness `ν` and p-norm `p`.
"""
function SeparableMaternCovarianceFunction(d::N where N<:Integer,λ::T where T<:Real,σ::T where T<:Real,ν::T where T<:Real,p::T where T<: Real)  
    d > 0 || throw(ArgumentError("dimension must be positive, got $(d)"))
    λ > 0 || throw(ArgumentError("correlation length λ of the Gaussian random field cannot be negative or zero!"))
    σ > 0 || throw(ArgumentError("marginal standard deviation σ of the Gaussian random field cannot be negative or zero!"))
    ν > 0 || throw(ArgumentError("smoothness ν of Gaussian random field cannot be negative or zero!"))
    p >= 1 || throw(ArgumentError("in p-norm, p must be greater than one!"))

    SeparableMaternCovarianceFunction{d,promote_type(typeof(λ),typeof(σ),typeof(ν),typeof(p))}(promote(λ,σ,ν,p)...)
end

"""
apply_covariance_function(M, x, y) where M<:AbstractMaternCovarianceFunction

Apply the Matern covariance function `M` to the points in the vectors `x` and `y`. In the one-dimensional case, `x` and `y` are vectors of real numbers. In the `d`-dimensional case, `x` and `y` are vectors of `d`-dimensional points.
"""
# TODO test d-dimensional implementation?
# TODO add example
function apply_covariance_function(K::M where M<:AbstractMaternCovarianceFunction, x::AbstractVector{T}, y::AbstractVector{T}) where T
    cov = zeros(eltype(T),length(x),length(y))
    for i in 1:length(x)
        for j in 1:length(y) # HACK we remove sigma^2 from the covariance matrix for numerical stability
            @inbounds cov[i,j] = x[i].==y[i] ? 1 : 2^(1-K.ν)/gamma(K.ν)*(sqrt(2*K.ν)*norm(x[i]-y[j],K.p)/K.λ).^K.ν.*besselk(K.ν,sqrt(2*K.ν)*norm(x[i]-y[j],K.p)/K.λ)
        end
    end

    return cov
end

# utilities
ndims(cov::AbstractMaternCovarianceFunction{d}) where d = d

# output formatting
function show(io::IO, M::MaternCovarianceFunction)
    print(io, "Mat\'ern covariance function in $(ndims(cov)) dimensions with correlation length λ = $(M.λ), marginal standard deviation σ = $(M.σ), smoothness ν = $(M.ν) and $(M.p)-norm")
end

function show(io::IO, M::SeparableMaternCovarianceFunction)
    print(io, "separable Mat\'ern covariance function in $(ndims(cov)) dimensions with correlation length λ = $(M.λ), marginal standard deviation σ = $(M.σ), smoothness ν = $(M.ν) and $(M.p)-norm")
end

## exponential covariance function ##
"""
ExponentialCovarianceFunction

Representation of an exponential covariance function.
"""
struct ExponentialCovarianceFunction{d,T} <: AbstractExponentialCovarianceFunction{d} where T<:Real
    λ::T
    σ::T
    p::T
end

"""
ExponentialCovarianceFunction(d, λ, σ, p)

Create an exponential covariance function in `d` dimensions with correlation length `λ`, marginal standard deviation `σ` and p-norm `p`.
"""
function ExponentialCovarianceFunction(d::N where N<:Integer,λ::T where T<:Real,σ::T where T<:Real,p::T where T<: Real)  
    d > 0 || throw(ArgumentError("dimension must be positive, got $(d)"))
    λ > 0 || throw(ArgumentError("correlation length λ of the Gaussian random field cannot be negative or zero!"))
    σ > 0 || throw(ArgumentError("marginal standard deviation σ of the Gaussian random field cannot be negative or zero!"))
    p >= 1 || throw(ArgumentError("in p-norm, p must be greater than one!"))

    ExponentialCovarianceFunction{d,promote_type(typeof(λ),typeof(σ),Float64,typeof(p))}(promote(λ,σ,p)...)
end

"""
SeparableExponentialCovarianceFunction

Representation of a separable  exponential covariance function.
"""
struct SeparableExponentialCovarianceFunction{d,p,T} <: AbstractExponentialCovarianceFunction{d} where T<:Real
    λ::T
    σ::T
end

"""
SeparableExponentialCovarianceFunction(d, λ, σ, p)

Create a separable exponential covariance function in `d` dimensions with correlation length `λ`, marginal standard deviation `σ` and p-norm `p`.
"""
function SeparableExponentialCovarianceFunction(d::N where N<:Integer,λ::T where T<:Real,σ::T where T<:Real,p::T where T<: Real)  
    d > 0 || throw(ArgumentError("dimension must be positive, got $(d)"))
    λ > 0 || throw(ArgumentError("correlation length λ of the Gaussian random field cannot be negative or zero!"))
    σ > 0 || throw(ArgumentError("marginal standard deviation σ of the Gaussian random field cannot be negative or zero!"))
    p >= 1 || throw(ArgumentError("in p-norm, p must be greater than one!"))

    SeparableExponentialCovarianceFunction{d,p,promote_type(typeof(λ),typeof(σ),Float64,typeof(p))}(promote(λ,σ)...)
end

"""
apply_covariance_function(M, x, y) where M<:AbstractExponentialCovarianceFunction

Apply the Expinential covariance function `M` to the points in the vectors `x` and `y`. In the one-dimensional case, `x` and `y` are vectors of real numbers. In the `d`-dimensional case, `x` and `y` are vectors of `d`-dimensional points.
"""
# TODO test d-dimensional implementation?
# TODO change this kernel!
# TODO add example
function apply_covariance_function(K::M where M<:AbstractExponentialCovarianceFunction, x::AbstractVector{T}, y::AbstractVector{T}) where T
    cov = zeros(eltype(T),length(x),length(y))
    for i in 1:length(x)
        for j in 1:length(y) # HACK we remove sigma^2 from the covariance matrix for numerical stability
            @inbounds cov[i,j] = x[i].==y[i] ? 1 : 2^(1-K.ν)/gamma(K.ν)*(sqrt(2*K.ν)*norm(x[i]-y[j],K.p)/K.λ).^K.ν.*besselk(K.ν,sqrt(2*K.ν)*norm(x[i]-y[j],K.p)/K.λ)
        end
    end

    return cov
end

# utilities
ndims(cov::AbstractExponentialCovarianceFunction{d}) where d = d
norm(cov::SeparableExponentialCovarianceFunction{d,p}) where {d,p} = p

# output formatting
function show(io::IO, E::ExponentialCovarianceFunction)
    print(io, "exponential covariance function in $(ndims(E)) dimensions with correlation length λ = $(M.λ), marginal standard deviation σ = $(M.σ) and $(M.p)-norm")
end

function show(io::IO, E::SeparableExponentialCovarianceFunction)
    print(io, "separable Mat\'ern covariance function in $(ndims(E)) dimension with correlation length λ = $(M.λ), marginal standard deviation σ = $(M.σ) and $(norm(M))-norm")
end

