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
    estimator.verbose && print_header(estimator, ϵ)

    # MC simulation is performed at Level(0)
    level = Level(0)
    push!(estimator, level)
    
    # obtain initial variance estimate
    if !contains_samples_at_index(estimator, level)
        sample!(estimator, level, nb_of_warm_up_samples(estimator))
    end

    # print status
    estimator.verbose && print_status(estimator)

    # evaluate optimal number of samples
    n = ceil.(Int,2/ϵ^2 * var(estimator, level))

    # print optimal number of samples
    estimator.verbose && print_optimal_nb_of_samples(estimator, n)

    # take additional samples
    n_due = n - nb_of_samples(estimator, level)
    n_due > 0 && sample!(estimator, level, n_due)

    # print convergence status
    estimator.verbose && print_convergence(estimator, true)
end

# bias
bias(estimator::Estimator{<:SL}) = 0.

# rates
α(estimator::Estimator{<:SL}) = NaN
β(estimator::Estimator{<:SL}) = NaN
γ(estimator::Estimator{<:SL}) = NaN
