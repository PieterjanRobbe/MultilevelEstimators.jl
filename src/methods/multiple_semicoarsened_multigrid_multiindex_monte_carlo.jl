## multiple_semicoarsened_multigrid_multiindex_monte_carlo.jl : Multiple Semi-Coarsened Multi-Index Monte Carlo method
#
# Implementation of the Multiple Semi-Coarsened Multi-Index Monte Carlo (MSG-MIMC) method.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## sample from exponential distribution with rate r
randexpr(r::Number, kwargs...) = randexp(kwargs...)/r

## varest ##
function varest(estimator::Estimator{<:MG, <:MC})
	idx = qoi_with_max_var(estimator)
	acc = accumulator(estimator, idx)
	length(acc) < 2 ? Inf : var(acc)/length(acc)
end

## set_weights ##
set_weights(estimator::Estimator{<:MG}) = set_weights(estimator, keys(estimator))
function set_weights(estimator::Estimator{<:MG}, index_set)
	for index in union(index_set, keys(estimator))
		set_weight(estimator, index, w(r(estimator), index))
	end
end

## r, p, w ##
function r(estimator::Estimator{<:MG})
	r = 1/2*(β(estimator) + γ(estimator))
	r = map(rᵢ->isnan(rᵢ) || rᵢ ≤ 0 ? 1.5 : rᵢ, r)
end
p(r, level::Level) = exp2(-r*level[1])
w(r, level::Level) = 1/p(r, level)
