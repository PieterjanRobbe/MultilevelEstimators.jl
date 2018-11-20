## monte_carlo.jl : standard Monte Carlo method
#
# Standard Monte Carlo method
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018
# TODO how to make comparison with QMC method?
#  - implement as qmc, with MersenneTwister and deticated var functions?
#  - deticated method?

## main routine ##
function _run(estimator::Estimator{<:SL, <:MC}, ϵ::T where {T<:Real})

    # print status
    verbose(estimator) && print_header(estimator, ϵ)

    # MC simulation is performed at Level(0)
    level = Level(0)
    push!(estimator, level)
    
    # obtain initial variance estimate
    if !contains_samples_at_index(estimator, level)
        sample!(estimator, level, nb_of_warm_up_samples(estimator))
    end

    # print status
    verbose(estimator) && print_status(estimator)

    # evaluate optimal number of samples
    n = ceil.(Int,2/ϵ^2 * var(estimator, level))

    # print optimal number of samples
    verbose(estimator) && print_optimal_nb_of_samples(estimator, n)

    # take additional samples
    n_due = n - nb_of_samples(estimator, level)
    n_due > 0 && sample!(estimator, level, n_due)

    # print convergence status
    verbose(estimator) && print_convergence(estimator, true)
end

# bias
bias(estimator::Estimator{<:SL}) = 0.0

# rates
α(estimator::Estimator{<:SL}) = NaN
β(estimator::Estimator{<:SL}) = NaN
γ(estimator::Estimator{<:SL}) = NaN

##########
cost(estimator::Estimator, index::Index) = total_work(estimator, index)/nb_of_samples(estimator, index)

qoi_with_max_var(estimator::Estimator{<:AbstractIndexSet, <:MC}) = argmax(map(q->sum(map(i->var(getindex(samples_diff(estimator, q), i)), keys(estimator))), 1:nb_of_qoi(estimator)))

var(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index) = var(samples_diff(estimator, qoi_with_max_var(estimator), index))
mean(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index) = mean(samples_diff(estimator, qoi_with_max_var(estimator), index))
varest(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index) = var(estimator, index)/nb_of_samples(estimator, index)
var(estimator::Estimator{<:AbstractIndexSet, <:MC}) = sum(var(estimator, index) for index in keys(estimator))
mean(estimator::Estimator{<:AbstractIndexSet, <:MC}) = sum(mean(estimator, index) for index in keys(estimator))
varest(estimator::Estimator{<:AbstractIndexSet, <:MC}) = sum(varest(estimator, index) for index in keys(estimator))

mse(estimator::Estimator) = varest(estimator) + bias(estimator)^2

rmse(estimator::Estimator) = sqrt(mse(estimator))
