## random_fields.jl : gaussian random fields covariance functions

## CovarianceFunction ##
"""
CovarianceFunction

Abstract covariance function type.
"""
abstract type CovarianceFunction end

## MaternCovarianceFunction ##
"""
MaternCovarianceFunction

Representation of a Matern covariance function.
"""
struct MaternCovarianceFunction{d,T} <: CovarianceFunction where {d,T<:Real}
    λ::T
    σ::T
    ν::T
    p::T
    is_separable::Bool
end

"""
MaternCovarianceFunction(d, λ, σ, ν, p, is_separable = true)

Create a Matern covariance function in d dimensions with correlation length `λ`, marginal standard deviation `σ`, smoothness `ν` and p-norm `p`.
"""
function MaternCovarianceFunction(d::N where N<:Integer,λ::T where T<:Real,σ::T where T<:Real,ν::T where T<:Real,p::T where T<: Real; is_separable::Bool = true)  
    d > 0 || throw(ArgumentError("dimension must be positive, got $(d)"))
    λ > 0 || throw(ArgumentError("correlation length λ of the Gaussian random field cannot be negative or zero!"))
    σ > 0 || throw(ArgumentError("marginal standard deviation σ of the Gaussian random field cannot be negative or zero!"))
    ν > 0 || throw(ArgumentError("smoothness ν of Gaussian random field cannot be negative or zero!"))
    p >= 1 || throw(ArgumentError("in p-norm, p must be greater than one!"))

    MaternCovarianceFunction{d,promote_type(typeof(λ),typeof(σ),typeof(ν),typeof(p))}(promote(λ,σ,ν,p)...,is_separable)
end

"""
apply_covariance_function(M, x, y) where M<:MaternCovarianceFunction

Apply the Matern covariance function `M` to the points in the vectors `x` and `y`. In the one-dimensional case, `x` and `y` are vectors of real numbers. In the `d`-dimensional case, `x` and `y` are vectors of `d`-dimensional points.
"""
# TODO test d-dimensional implementation?
function apply_covariance_function(K::MaternCovarianceFunction,x::AbstractVector{T},y::AbstractVector{T}) where T
    cov = zeros(eltype(T),length(x),length(y))
    for i in 1:length(x)
        for j in 1:length(y) # HACK we remove sigma^2 from the covariance matrix for numerical stability
            @inbounds cov[i,j] = x[i].==y[i] ? 1 : 2^(1-K.ν)/gamma(K.ν)*(sqrt(2*K.ν)*norm(x[i]-y[j],K.p)/K.λ).^K.ν.*besselk(K.ν,sqrt(2*K.ν)*norm(x[i]-y[j],K.p)/K.λ)
        end
    end

    return cov
end

# output formatting
function show(io::IO, M::MaternCovarianceFunction)
    print(io, "Matern kernel with correlation length λ = $(M.λ), marginal standard deviation σ = $(M.σ), smoothness ν = $(M.ν) and $(M.p)-norm")
end

## ExponentialCovarianceFunction ##
"""
ExponentialCovarianceFunction(d, λ, σ, p, is_separable = true)

Create an exponential covariance function in d dimensions with correlation length `λ`, marginal standard deviation `σ` and p-norm `p`.
"""
function ExponentialCovarianceFunction(d::N where N<:Integer,λ::T where T<:Real,σ::T where T<:Real,p::T where T<: Real; is_separable::Bool = true)  
    d > 0 || throw(ArgumentError("dimension must be positive, got $(d)"))
    λ > 0 || throw(ArgumentError("correlation length λ of the Gaussian random field cannot be negative or zero!"))
    σ > 0 || throw(ArgumentError("marginal standard deviation σ of the Gaussian random field cannot be negative or zero!"))
    p >= 1 || throw(ArgumentError("in p-norm, p must be greater than one!"))

    MaternCovarianceFunction{d,promote_type(typeof(λ),typeof(σ),Float64,typeof(p))}(promote(λ,σ,0.5,p)...,is_separable)
end

