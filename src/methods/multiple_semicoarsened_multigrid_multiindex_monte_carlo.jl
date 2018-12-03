## multiple_semicoarsened_multigrid_multiindex_monte_carlo.jl : Multiple Semi-Coarsened Multi-Index Monte Carlo method
#
# Implementation of the Multiple Semi-Coarsened Multi-Index Monte Carlo (MSG-MIMC) method.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## sample from exponential distribution with rate r
randexpr(r::Number, kwargs...) = randexp(kwargs...)/r

## varest ##
function varest(estimator::Estimator{I, <:MC}) where I<:MG
    idx = qoi_with_max_var(estimator)
    @show N = length(samples(estimator, idx, zero(eltype(I))))
    #@show N = length(samples(estimator, idx, one(eltype(I))))
    if N == 1
        return Inf
    else
        Ys = zeros(N)
        for index in keys(estimator)
            @show w = weight_factor(estimator, index)
            ns = length(samples(estimator, idx, index))
            Ys[N-ns+1:N] .+= w*samples_diff(estimator, idx, index)
        end
        return var(Ys)/N
    end
end

## compute weight factor
function weight_factor(estimator::Estimator, index::Index)
    r = 1/2*(β(estimator) + γ(estimator))
    r = map(rᵢ->isnan(rᵢ) || rᵢ ≤ 0 ? 1.5 : rᵢ, r)
    w(r, index)
end

p(r, level::Level) = exp2(-r*level[1])
w(r, level::Level) = 1/p(r, level)
