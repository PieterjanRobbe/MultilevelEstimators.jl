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
        println("plain Var0 is $(var(estimator, zero(eltype(I)))/length(samples_diff(estimator, idx, zero(eltype(I)))))")
        return var(Ys)/N
    end
end

## compute weight factor
function weight_factor(estimator::Estimator, index::Index)
    r = 1/2*(β(estimator) + γ(estimator))
    r = map(rᵢ->isnan(rᵢ) || rᵢ ≤ 0 ? 1.5 : rᵢ, r)
    #r = log(2)*1.5
    @show r
    #prod(map(i->exp(r[i]*index[i])*(1-exp(-r[i]))/r[i], 1:ndims(estimator)))
    #prod(map(i->exp(r[i]*(index[i])), 1:ndims(estimator)))
        p = Dict(index=>2.0^(-r*index[1]) for index in keys(estimator))
        display(p)
        println("")
        sum(values(p))/p[index]
end
