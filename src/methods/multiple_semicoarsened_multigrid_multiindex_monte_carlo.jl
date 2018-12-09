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
    # compute "optimal" r
    r = isempty(keys(estimator)) ? ntuple(i->1.5, ndims(estimator)) : 1/2*(β(estimator) + γ(estimator))
    r = map(rᵢ->isnan(rᵢ) || rᵢ≤0 ? 1.5 : rᵢ, r)
    # limit r to avoid catastrophic failure
    for index in keys(estimator)
        MQ = mean(estimator, ntuple(i->0, ndims(estimator)))
        MY = mean(estimator, index)
        max_r = log(abs(MQ/MY))./index
        r = map(i->max_r[i] ≤ 0 || !isfinite(max_r[i]) ? r[i] : min(max_r[i], r[i]), 1:ndims(estimator))
    end
    r
end
p(r, index::Index) = prod(exp2.(.-r.*index))
w(r, index::Index) = inv.(p(r, index))
