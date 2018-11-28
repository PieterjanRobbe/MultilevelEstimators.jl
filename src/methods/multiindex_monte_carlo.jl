## multiindex_monte_carlo.jl : Multi-Index Monte Carlo method
#
# Implementation of the Multi-Index Monte Carlo (MIMC) method.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## rates ##
for (f, g, sgn) in zip((:α, :β, :γ), (:mean, :var, :cost), (-1, -1, 1))
	eval(
		 quote
			 $f(estimator::Estimator{<:AbstractMI}) = $sgn.*getindex.(broadcast(i->$(Symbol("rates_", f))(estimator, i), 1:ndims(estimator)), 2)
			 $(Symbol("rates_", f))(estimator::Estimator{<:AbstractMI}, dir::Integer) = $(Symbol("rates_", f))(estimator, (maximum(getindex.(keys(estimator), dir)) + 1) .* ntuple(i->i==dir, ndims(estimator)), dir)
			 function $(Symbol("rates_", f))(estimator::Estimator{<:AbstractMI}, idx::Index, dir::Integer)
				 m = idx[dir] - 1
				 if m < 2
					 return (NaN, NaN)
				 else
					 x = m:-1:1
					 y = map(xᵢ -> $g(estimator, idx .- xᵢ.*ntuple(i->i==dir, ndims(estimator))), 1:m)
					 idcs = .!isnan.(y)
					 θ = interp1(view(x, idcs), log2.(abs.(view(y, idcs))))
					 return tuple(θ...)
				 end
			 end
		 end)
end

## regression ##
function interpd(estimator::Estimator{I}, f::Function) where I<:AbstractMI
	idx_set = all_keys(estimator)
	A = [getindex(idx_set[i], j) for i in 2:length(idx_set), j in 1:ndims(estimator)]
	A = hcat(ones(eltype(A), size(A, 1)), A)
	y = map(index->log2(f(estimator, index)), view(idx_set, 2:length(idx_set)))
	A\y
end

for (f, sym) in zip([:var, :cost], [:β, :γ])
	eval(
		 quote
			 function $(Symbol("_regress_", f))(estimator::Estimator{<:MI, <:MC}, index::Index)
				 p = broadcast(dir->$(Symbol("rates_", sym))(estimator, index, dir), 1:ndims(estimator))
				 estimates = broadcast(dir->2^(p[dir][1]+index[dir]*p[dir][2]), 1:ndims(estimator))
				 estimate = mean(filter(!isnan, estimates))
				 if isnan(estimate)
					 p = interpd(estimator, $f)
					 return 2 .^(p[1]+sum(p[2:end].*index))
				 else
				   	 return estimate
				 end
			 end
		 end)
end

regress_var(estimator::Estimator{<:MI, <:MC}, index::Index) = _regress_var(estimator::Estimator{<:MI, <:MC}, index::Index)
regress_cost(estimator::Estimator{<:MI, <:MC}, index::Index) = cost_model(estimator) isa EmptyFunction ? _regress_cost(estimator::Estimator{<:MI, <:MC}, index::Index) : cost_model(estimator, index)

function _regress_nb_of_samples(estimator::Estimator{<:AbstractIndexSet, <:MC}, index_set::AbstractVector{<:Index}, ϵ::Real, θ::Real)
	vars = Dict(index=>regress_var(estimator, index) for index in index_set)
	costs = Dict(index=>regress_cost(estimator, index) for index in index_set)
	Σ_estimate = Σ(estimator)
	for index in keys(vars)
		Σ_estimate += sqrt(vars[index] * costs[index])
	end
	Dict(index=>max(2, min(optimal_nb_of_samples(ϵ, θ, vars[index], costs[index], Σ_estimate), nb_of_warm_up_samples(estimator))) for index in index_set)
end

## bias ##
bias(estimator::Estimator{<:AbstractMI}) = bias(estimator, sz(estimator))
function bias(estimator::Estimator{<:AbstractMI}, sz::Integer)
	if !isempty(boundary(estimator, sz+1) ∩ keys(estimator)) && !robustify_bias_estimate(estimator)
		return abs(sum(broadcast(i->mean(estimator, i), boundary(estimator, sz+1))))
	else
		x = 1:sz
		y = [log2(abs(sum(broadcast(i->mean(estimator, i), boundary(estimator, xᵢ))))) for xᵢ in x]
		p = interp1(x, y)
		return 2^(p[1]+(sz+1)*p[2])
	end
end
