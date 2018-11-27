## multilevel_monte_carlo.jl : Multilevel Monte Carlo method
#
# Implementation of the Multilevel Monte Carlo (MLMC) method.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## converged ##
converged(estimator::Estimator, ϵ::Real, θ::Real) = ( bias(estimator)^2 ≤ (1-θ)*ϵ^2 || mse(estimator) ≤ ϵ^2 )

## MSE and RMSE ##
mse(estimator::Estimator) = varest(estimator) + bias(estimator)^2
rmse(estimator::Estimator) = sqrt(mse(estimator))

## qoi_with_max_var ##
qoi_with_max_var(estimator::Estimator{<:AbstractIndexSet, <:MC}) = argmax(map(q->sum(map(i->var(getindex(samples_diff(estimator, q), i)), keys(estimator))), 1:nb_of_qoi(estimator)))

## mean0 and var0 at index ##
for f in [:mean, :var]
	eval(
		 quote
			 $(Symbol(f, 0))(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index) = $f(samples(estimator, qoi_with_max_var(estimator), index))
		 end)
end

## cost at index ##
cost(estimator::Estimator, index::Index) = total_work(estimator, index)/nb_of_samples(estimator, index)

## mean and var at index ##
for f in [:mean, :var]
	eval(
		 quote
			 $f(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index) = $f(samples_diff(estimator, qoi_with_max_var(estimator), index))
		 end)
end

## varest at index ##
varest(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index) = var(estimator, index)/nb_of_samples(estimator, index)

## mean, var, varest ##
for f in [:mean, :var, :varest]
	eval(
		 quote
			 $f(estimator::Estimator{<:AbstractIndexSet, <:MC}) = sum($f(estimator, index) for index in keys(estimator))
		 end)
end

## rates ##
function interp1(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
	A = hcat(ones(eltype(x), length(x)), x)
	A\y
end

for (f, g, sgn) in zip((:α, :β, :γ), (:mean, :var, :cost), (-1, -1, 1))
	eval(
		 quote
			 $f(estimator::Estimator{<:AbstractML}) = $sgn*$(Symbol("rates_", f))(estimator)[2]
			 $(Symbol("rates_", f))(estimator::Estimator{<:AbstractML}) = $(Symbol("rates_", f))(estimator, sz(estimator))
			 $(Symbol("rates_", f))(estimator::Estimator{<:AbstractML}, max_level::Integer) = $(Symbol("rates_", f))(estimator, max_level-1, max_level)
			 function $(Symbol("rates_", f))(estimator::Estimator{<:AbstractML}, start_level::Integer, max_level::Integer)
				 if max_level < 2
					 return (NaN, NaN)
				 else
					 x = 1:max_level
					 y = map(xᵢ -> $g(estimator, Level(xᵢ)), x)
					 idcs = .!isnan.(y)
					 θ = interp1(view(x, idcs), log2.(abs.(view(y, idcs))))
					 return tuple(θ...)
				 end
			 end
		 end)
end

## regression ##
function regress_nb_of_samples(estimator::Estimator{<:ML, <:MC}, level::Level, ϵ::Real, θ::Real)
	p1 = rates_β(estimator)
	var_estimate = 2^(p1[1]+level[1]*p1[2])
	p2 = rates_γ(estimator)
	cost_estimate = 2^(p2[1]+level[1]*p2[2])
	Σ_estimate = Σ(estimator)
	Σ_estimate += sqrt(var_estimate * cost_estimate)
	max(2, min(optimal_nb_of_samples(ϵ, θ, var_estimate, cost_estimate, Σ_estimate), nb_of_warm_up_samples(estimator)))
end

## bias ##
bias(estimator::Estimator{<:AbstractML}) = bias(estimator, sz(estimator))
function bias(estimator::Estimator{<:AbstractML}, sz::Integer)
	start_level = robustify_bias_estimate(estimator) ? 1 : sz - 1 
	p = rates_α(estimator, start_level, sz)
	2^(p[1]+(sz+1)*p[2])
end

## compute optimal value of MSE splitting parameter ##
function compute_splitting(estimator::Estimator, ϵ::Real)
	bias_estimate = bias(estimator, max_sz(estimator))
	θ = splitting(estimator)
	isnan(bias_estimate) ? θ : min(0.99, max(θ, 1-bias_estimate^2/ϵ^2))
end

## optimal nb of samples ##
Σ(estimator::Estimator) = sum(sqrt.(map(index -> var(estimator, index) * cost(estimator, index), keys(estimator))))
optimal_nb_of_samples(estimator::Estimator, index::Index, ϵ::Real, θ::Real) = optimal_nb_of_samples(ϵ, θ, var(estimator, index), cost(estimator, index), Σ(estimator))
optimal_nb_of_samples(ϵ::Real, θ::Real, var_estimate::Real, cost_estimate::Real, Σ_estimate::Real) = ceil(Int, 1/(θ*ϵ^2) * sqrt(var_estimate/cost_estimate) * Σ_estimate)

## new index set ##
new_index_set(estimator::Estimator{<:Union{SL, AbstractML}}, n::Integer) = Set((Level(n),))
