## NonIsotropicMatern.jl : stationary non-isotropic Matern covariance function
#
# Returns a stationary non-isotropic Matern covariance function.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

struct NonIsotropicMatern{T} <: AnisotropicCovarianceStructure{T}
	λ₁::T
	λ₂::T
	σ::T
	ν::T
end

NonIsotropicMatern(λ₁::Real, λ₂::Real, ν::Real) = NonIsotropicMatern(λ₁, λ₂, one(typeof(λ₁)), ν)

rnorm(n::NonIsotropicMatern, x::Vector{<:Real}) = sqrt(x[1]^2/n.λ₁^2+x[2]^2/n.λ₂^2)

function GaussianRandomFields.apply(n::NonIsotropicMatern, x::Vector{<:Real})
	r = rnorm(n, x)
	if iszero(r)
		one(eltype(x))
	else
		tmp = sqrt(2 * n.ν) * r
		2^(1 - n.ν) / gamma(n.ν) * tmp^n.ν * besselk(n.ν, tmp)
	end
end

GaussianRandomFields.shortname(::NonIsotropicMatern) = "Non-isotropic Mat\u00E9rn"
