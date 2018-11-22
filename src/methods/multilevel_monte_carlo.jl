## multilevel_monte_carlo.jl : Multilevel Monte Carlo method
#
# Implementation of the Multilevel Monte Carlo method.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## main routine ##
function _run(estimator::Estimator{<:ML, <:MC}, ϵ::T where {T<:Real})

    # print status
	verbose(estimator) && print_header(estimator, ϵ)

    # start level is 0
    level = Level(0)

    # MSE splitting parameter
	θ = splitting(estimator)

    # main loop
	while !converged(estimator, level, ϵ, θ) 

        # obtain initial variance estimate
		if !contains_samples_at_index(estimator, level)
			if do_regression(estimator) &&  level > 2
				n = regress_nb_of_samples(estimator, level, ϵ, θ)
			else
				n = nb_of_warm_up_samples(estimator)
			end
			sample!(estimator, level, n)
		end

        # add new level to the index set
        push!(estimator, level)

        # print status
		verbose(estimator) && print_status(estimator)

        # value of the MSE splitting parameter
		θ = do_mse_splitting(estimator) ? compute_splitting(estimator, ϵ) : splitting(estimator)

        # evaluate optimal number of samples
		n_opt = Dict(τ => optimal_nb_of_samples(estimator, τ, ϵ, θ) for τ in keys(estimator))

        # print optimal number of samples
		verbose(estimator) && print_nb_of_samples(estimator, n_opt)

        # take additional samples
        for τ in keys(estimator)
			n_due = n_opt[τ] - nb_of_samples(estimator, τ)
            n_due > 0 && sample!(estimator, τ, n_due)
        end

        # show status
		verbose(estimator) && print_status(estimator) && print_mse_analysis(estimator, ϵ, θ)

        # increase level
        level = level + 1 

        # check if the new level exceeds the maximum level
		if !converged(estimator, level, ϵ, θ) && ( level > max_index_set_param(estimator) ) 
			verbose(estimator) && warn_max_level(estimator)
            break
        end
    end

    # print convergence status
    estimator.verbose && print_convergence(estimator,converged)
end

## converged ##
converged(estimator::Estimator, level::Level, ϵ::Real, θ::Real) = ( level ≥ 2 ) && ( bias(estimator)^2 ≤ (1-θ)*ϵ^2 || mse(estimator) ≤ ϵ^2 )

## cost ##
cost(estimator::Estimator, index::Index) = total_work(estimator, index)/nb_of_samples(estimator, index)

## qoi_with_max_var ##
qoi_with_max_var(estimator::Estimator{<:AbstractIndexSet, <:MC}) = argmax(map(q->sum(map(i->var(getindex(samples_diff(estimator, q), i)), keys(estimator))), 1:nb_of_qoi(estimator)))

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

## MSE and RMSE ##
mse(estimator::Estimator) = varest(estimator) + bias(estimator)^2
rmse(estimator::Estimator) = sqrt(mse(estimator))

## rates ##
function interp1(x::Vector{<:Real}, y::Vector{<:Real})
	A = hcat(ones(eltype(x), length(x)), x)
    A\y
end

for (f, g, sgn) in zip((:α, :β, :γ), (:mean, :var, :cost), (-1, -1, 1))
	eval(
		 quote
			 $f(estimator::Estimator{<:AbstractML}) = $sgn*$(Symbol("rates_", f))(estimator)[2]
			 $(Symbol("rates_", f))(estimator::Estimator{<:AbstractML}) = $sgn*$(Symbol("rates_", f))(estimator, first(maximum(keys(estimator))))[2]
			 $(Symbol("rates_", f))(estimator::Estimator{<:AbstractML}, max_level::Integer) = $(Symbol("rates_", f))(estimator, max_level-1, max_level)
			 function $(Symbol("rates_", f))(estimator::Estimator{<:AbstractML}, start_level::Integer, max_level::Integer)
				 if max_level < 2
					 return (NaN, NaN)
				 else
					 x = 1:max_level
					 y = map(xᵢ -> $g(estimator, Level(xᵢ)) for xᵢ in x)
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
	max(2, min(optimal_nb_of_samples(estimator, level, ϵ, θ), nb_of_warm_up_samples(estimator)))
end

## bias ##
bias(estimator::Estimator{<:AbstractML}) = bias(estimator, maximum(keys(estimator)))
function bias(estimator::Estimator{<:AbstractML}, max_level::Integer)
	start_level = robustify_bias_estimate(estimator) ? 1 : max_level-1 
	p = rates_α(estimator, start_level, max_level)
    2^(p[1]+(max_level+1)*p[2])
end

## compute optimal value of MSE splitting parameter ##
function compute_splitting(estimator::Estimator, ϵ::Real)
	bias_estimate = bias(estimator, max_level_where_samples_are_taken(estimator)[1])
	θ = splitting(estimator)
	isnan(bias_estimate) ? θ : min(0.99, max(θ, 1-bias_estimate^2/ϵ^2))
end

## optimal nb of samples ##
Σ(estimator::Estimator) = sum(sqrt.(map(index -> var(estimator, index) * cost(estimator, index), keys(estimator))))
optimal_nb_of_samples(estimator::Estimator, index::Index, ϵ::Real, θ::Real) = ceil(Int, 1/(θ*ϵ^2) * sqrt(var(estimator, index) / cost(estimator, index)) * Σ(estimator))
