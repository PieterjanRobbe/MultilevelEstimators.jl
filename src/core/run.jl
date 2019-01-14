## run.jl : entry point method for simulating an Estimator
#
# Runs the Estimator for the given tolerance, and returns a History object.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

"""
    run(estimator, ϵ)

Run the estimator and compute the expected value of the quantity of interest up to the given tolerance(s) `ϵ`.

# Examples
```jldoctest
julia>

```
"""
function run(estimator::Estimator, tols::AbstractVector{<:Real})

    # input checking
    check_larger_than("run", tols, 0, "supplied tolerance(s)")

    # make history
    history = History()

    # run the sequence of tolerances
    for tol in tols
        _run(estimator,tol)
        push!(history, estimator, tol) # log the results in history
        clear(estimator) # prepare new run
    end

    return history
end

run(estimator::Estimator, tol::Real) = run(estimator, get_tols(estimator, tol))

## Main routine ##
function _run(estimator::Estimator{T, <:MC}, ϵ::Real) where T<:AbstractIndexSet

    # print status
    verbose(estimator) && print_header(estimator, ϵ)

    # start "level" is 0
    L = 0

    # initial MSE splitting parameter
    θ = min_splitting(estimator)

    # main loop
    is_converged = false
    while !is_converged

        # print level
        verbose(estimator) && print_level(estimator, L)

        # update index set
        index_set = boundary(estimator, L)
        verbose(estimator) && print_index_set(estimator, index_set)

        # obtain initial variance estimate
        ns = regress_nb_of_samples(estimator, index_set, ϵ, θ, L) 
        for index in index_set
            if !contains_samples_at_index(estimator, index)
                sample!(estimator, index, ns[index])
            end
        end

        # add new indices to the index set
        for index in index_set
            push!(estimator, index)
        end
        set_sz(estimator, L)

        # print status
        verbose(estimator) && print_status(estimator)

        # update value of the MSE splitting parameter
        θ = do_mse_splitting(estimator) ? compute_splitting(estimator, ϵ) : min_splitting(estimator)

        # evaluate optimal number of samples
        n_opt = Dict(τ => optimal_nb_of_samples(estimator, τ, ϵ, θ) for τ in keys(estimator))

        # print optimal number of samples
        verbose(estimator) && print_optimal_nb_of_samples(estimator, n_opt)

        # take additional samples if required
        for τ in keys(estimator)
            n_due = n_opt[τ] - nb_of_samples(estimator, τ)
            n_due > 0 && sample!(estimator, τ, n_due)
        end
 
        # check for convergence
        is_converged = L > 2 && converged(estimator, ϵ, θ)
        if T <: SL
            is_converged = true
        else
            # show status
            verbose(estimator) && print_status(estimator)
            verbose(estimator) && L ≥ 2 && print_mse_analysis(estimator, ϵ, θ)

            # check if the new level exceeds the maximum level
            if !is_converged && ( sz(estimator) ≥ max_index_set_param(estimator) ) 
                verbose(estimator) && warn_max_level(estimator)
                is_converged = true
            else
                # update level
                L += 1
            end
        end
    end

    # print convergence status
    verbose(estimator) && print_convergence(estimator, converged(estimator, ϵ, θ))
end
