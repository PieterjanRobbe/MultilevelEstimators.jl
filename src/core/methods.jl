## methods.jl : refactored methods from run.jl
#
# Implements refactored methods from run.jl.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2019

ndims(::Estimator{<:AbstractIndexSet{d}}) where d = d 

get_index_set(estimator::Estimator, sz) = get_index_set(estimator.index_set, sz)

get_tols(estimator::Estimator, tol::T) where T<:Real = estimator[:continuate] ? estimator[:continuation_mul_factor].^(estimator[:nb_of_tols]-1:-1:0)*tol : T[tol] 

mse(estimator::Estimator) = varest(estimator) + bias(estimator)^2

rmse(estimator::Estimator) = sqrt(mse(estimator))

converged(estimator::Estimator, ϵ::Real, θ::Real) = ( bias(estimator)^2 ≤ (1-θ)*ϵ^2 || mse(estimator) ≤ ϵ^2 )

#
# inspector functions: mean, var, varest...
#
qoi_with_max_var(estimator::Estimator{<:AbstractIndexSet, <:MC}) = argmax(map(q->sum(map(i->var(getindex(samples_diff(estimator, q), i + one(i))), keys(estimator))), 1:estimator[:nb_of_qoi]))

cost(estimator::Estimator, index::Index) = estimator[:cost_model] isa EmptyFunction ? time(estimator, index) : work(estimator, index)

for f in [:mean, :var]
	@eval begin
		$f(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index) = $f(samples_diff(estimator, qoi_with_max_var(estimator), index))
		$(Symbol(f, 0))(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index) = $f(samples(estimator, qoi_with_max_var(estimator), index))
		$f(estimator::Estimator{<:AbstractIndexSet, <:MC}) = sum($f(estimator, index) for index in keys(estimator))
	end
end

varest(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index) = var(estimator, index)/nb_of_samples(estimator, index)

varest(estimator::Estimator{<:AbstractIndexSet, <:MC}) =  sum(varest(estimator, index) for index in keys(estimator))

#
# rates 
#
for (f, g, sgn) in zip([:α, :β, :γ], [:mean, :var, :cost], [-1, -1, 1])
	@eval begin 

		$f(estimator::Estimator{<:SL}) = nothing
		
		$f(estimator::Estimator{<:AbstractIndexSet}) = $sgn.*getindex.(broadcast(i->$(Symbol("rates_", f))(estimator, i), 1:ndims(estimator)), 2)

		$(Symbol("rates_", f))(estimator::Estimator{<:AbstractML}) = $(Symbol("rates_", f))(estimator, 1) 
		
		$(Symbol("rates_", f))(estimator::Estimator{<:AbstractIndexSet}, dir::Integer) = $(Symbol("rates_", f))(estimator, (maximum(getindex.(keys(estimator), dir)) + 1) * Index(ntuple(i -> i == dir, ndims(estimator))), dir)
		
		function $(Symbol("rates_", f))(estimator::Estimator{<:AbstractIndexSet}, idx::Index, dir::Integer)
			m = idx[dir] - 1
			if m < 2
				return (NaN, NaN)
			else
				x = m:-1:1
				y = map(xᵢ -> $g(estimator, idx - xᵢ * Index(ntuple(i -> i == dir, ndims(estimator)))), 1:m)
				idcs = .!isnan.(y)
				θ = interp1(view(x, idcs), log2.(abs.(view(y, idcs))))
				return tuple(θ...)
			end
		end
		
		function $(Symbol("_regress_", g))(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index)
			p = broadcast(dir->$(Symbol("rates_", f))(estimator, index, dir), 1:ndims(estimator))
			estimates = broadcast(dir->2^(p[dir][1]+index[dir]*p[dir][2]), 1:ndims(estimator))
			estimate = mean(filter(!isnan, estimates))
			if isnan(estimate)
				p = interp($g, estimator)
				return 2 .^(p[1]+sum(p[2:end].*index.I))
			else
				return estimate
			end
		end
	end
end

#
# regression
#
function interp1(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
	A = hcat(ones(eltype(x), length(x)), x)
	A\y
end

function interp(f::Function, estimator::Estimator{I}) where I<:AbstractMI
	idx_set = filter(i -> !isempty(samples(estimator)[1][i]), CartesianIndices(size(samples(estimator)[1])))
	A = [i == 0 ? 1 : getindex(index - one(index), i) for index in idx_set, i in 0:ndims(estimator)]	
	y = map(i -> log2(f(estimator, i - one(i))), idx_set)
	A\y
end

regress_mean(estimator, index) = _regress_mean(estimator, index)
regress_var(estimator, index) = _regress_var(estimator, index)
regress_cost(estimator, index) = estimator[:cost_model] isa EmptyFunction ? _regress_cost(estimator, index) : estimator[:cost_model](index)

function regress_nb_of_samples(estimator::Estimator, index_set::AbstractVector{<:Index}, ϵ::Real, θ::Real, L::Integer)
	if estimator[:do_regression] && L > 2
		return _regress_nb_of_samples(estimator, index_set, ϵ, θ)
	else
		return Dict(index => estimator[:nb_of_warm_up_samples] for index in index_set)
	end
end

# TODO define regress_var en regress_cost
function _regress_nb_of_samples(estimator::Estimator{<:AbstractIndexSet, <:MC}, index_set::AbstractVector{<:Index}, ϵ::Real, θ::Real)
	vars = Dict(index => regress_var(estimator, index) for index in index_set)
	costs = Dict(index => regress_cost(estimator, index) for index in index_set)
	Σ_estimate = Σ(estimator)
	for index in keys(vars)
		Σ_estimate += sqrt(vars[index] * costs[index])
	end
	Dict(index => begin
			 n_opt = optimal_nb_of_samples(ϵ, θ, vars[index], costs[index], Σ_estimate)
			 max(2, min(n_opt, estimator[:nb_of_warm_up_samples])) 
		 end for index in index_set)
end

function compute_splitting(estimator::Estimator, ϵ::Real)
	bias_estimate = bias(estimator, max_sz(estimator))
	a = estimator[:min_splitting]
	b = estimator[:max_splitting]

	splitting = 1 - bias_estimate^2/ϵ^2
	isnan(splitting) ? a : min(b, max(a, splitting))
end

Σ(estimator::Estimator) = sum(sqrt.(map(index -> var(estimator, index) * cost(estimator, index), keys(estimator))))

optimal_nb_of_samples(estimator::Estimator, index::Index, ϵ::Real, θ::Real) = optimal_nb_of_samples(ϵ, θ, var(estimator, index), cost(estimator, index), Σ(estimator))

optimal_nb_of_samples(ϵ::Real, θ::Real, var_estimate::Real, cost_estimate::Real, Σ_estimate::Real) = ceil(Int, 1/(θ*ϵ^2) * sqrt(var_estimate/cost_estimate) * Σ_estimate)

#
# bias computation
#
boundary(estimator::Estimator, cntr) = setdiff(get_index_set(estimator, cntr), get_index_set(estimator, cntr-1))

bias(estimator::Estimator{<:SL}) = 0.0

bias(estimator::Estimator) = bias(estimator, sz(estimator))

function bias(estimator::Estimator, sz::Integer)
	if !isempty(boundary(estimator, sz+1) ∩ keys(estimator)) && !robustify_bias_estimate(estimator)
		return abs(sum(broadcast(i -> mean(estimator, i), boundary(estimator, sz+1))))
	else
		x = 1:sz
		y = Float64[log2(abs(sum(broadcast(i -> mean(estimator, i), boundary(estimator, xᵢ))))) for xᵢ in x]
		p = interp1(x, y)
		return 2^(p[1]+(sz+1)*p[2])
	end
end
