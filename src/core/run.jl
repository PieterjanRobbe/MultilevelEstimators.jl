## run.jl : entry point method for simulating an Estimator
#
# Runs the Estimator for the given tolerance, and returns a History object.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2019

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
    check_larger_than("run", tols, "supplied tolerance(s)", 0)

    # make history
    history = History()

    # run the sequence of tolerances
    for tol in tols
        clear(estimator) # prepare new run
        t = @elapsed _run(estimator,tol)
        push!(history, estimator, tol, t) # log the results in history
    end

    return history
end

run(estimator::Estimator, tol::Real) = run(estimator, get_tols(estimator, tol))

## Main routine ##
function _run(estimator::Estimator{T, <:MC}, ϵ::Real) where T<:AbstractIndexSet

    # print status
	estimator[:verbose] && print_header(estimator, ϵ)

    # start "level" is 0
    L = 0

    # initial MSE splitting parameter
	θ = estimator[:min_splitting]

    # main loop
    is_converged = false
    while !is_converged

        # print level
        estimator[:verbose] && print_level(estimator, L)

        # update index set
        index_set = boundary(estimator, L)
        estimator[:verbose] && print_index_set(estimator, index_set)

        # obtain initial variance estimate
        ns = regress_nb_of_samples(estimator, index_set, ϵ, θ, L) 
        for index in index_set
            if !has_samples_at_index(estimator, index)
                sample!(estimator, index, ns[index])
            end
        end

        # add new indices to the index set
        for index in index_set
            push!(estimator, index)
        end
        set_sz(estimator, L)

        # print status
        estimator[:verbose] && print_status(estimator)

        # update value of the MSE splitting parameter
		θ = estimator[:do_mse_splitting] ? compute_splitting(estimator, ϵ) : estimator[:min_splitting]

        # evaluate optimal number of samples
        n_opt = Dict(τ => optimal_nb_of_samples(estimator, τ, ϵ, θ) for τ in keys(estimator))

        # print optimal number of samples
        estimator[:verbose] && print_optimal_nb_of_samples(estimator, n_opt)

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
            estimator[:verbose] && print_status(estimator)
            estimator[:verbose] && L ≥ 2 && print_mse_analysis(estimator, ϵ, θ)

            # check if the new level exceeds the maximum level
			if !is_converged && ( sz(estimator) ≥ estimator[:max_index_set_param] ) 
                estimator[:verbose] && warn_max_level(estimator)
                is_converged = true
            else
                # update level
                L += 1
            end
        end
    end

    # print convergence status
    estimator[:verbose] && print_convergence(estimator, converged(estimator, ϵ, θ))
end
