## run.jl : entry point method for simulating an Estimator
#
# Runs the Estimator for the given tolerance, and returns a History object.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for
# Multilevel Monte Carlo Methods (c) Pieterjan Robbe, 2019

"""
run(estimator::Estimator, ε::Real)
run(estimator::Estimator, ε::Vector{<:Real})

Run the estimator and compute the expected value of the quantity of interest up to the given tolerance(s) `ε`. Returns a [`History`](@ref)-object that contains usefull diagnostics about the simulation. 

# Examples

An example detailing how to use MultilevelEstimators is provided in the [example in the documentation](@ref Example).

See also: [`Estimator`](@ref), [`History`](@ref)
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
        update_history!(history, estimator, tol, t) # log the results in history
    end

    return history
end

run(estimator::Estimator, tol::Real) = run(estimator, get_tols(estimator, tol))

## Main routine ##
function _run(estimator::Estimator{T1, T2}, ϵ::Real) where {T1<:AbstractIndexSet, T2<:AbstractSampleMethod}

    # print status
    estimator[:verbose] && print_header(estimator, ϵ)

    # start "level" is 0
    L = 0

    # initial MSE splitting parameter
    θ = T1 <: SL ? 0.5 : estimator[:min_splitting]

    # restart flag
    flag = false

    # main loop
    is_converged = false
    while !is_converged

        # print level
        estimator[:verbose] && print_level(estimator, L)

        # update index set
        index_set = new_index_set(estimator, L)
        estimator[:verbose] && print_index_set(estimator, index_set)

        # obtain initial variance estimate
        ns = regress_nb_of_samples(estimator, index_set, ϵ, θ, L) 
        for index in index_set
            if !has_samples_at_index(estimator, index)
                flag = sample!(estimator, index, ns[index])
                flag && break
            end
        end

        flag && break

        # add new indices to the index set
        for index in index_set
            push!(estimator, index)
        end
        set_sz(estimator, L)

        # print status
        estimator[:verbose] && print_status(estimator)

        # update value of the MSE splitting parameter
        θ = T1 <: SL ? 0.5 : estimator[:do_mse_splitting] ? compute_splitting(estimator, ϵ) : estimator[:min_splitting]

        if T2 <: MC

            # evaluate optimal number of samples
            n_opt = Dict(τ => optimal_nb_of_samples(estimator, τ, ϵ, θ) for τ in keys(estimator))

            # print optimal number of samples
            estimator[:verbose] && print_optimal_nb_of_samples(estimator, n_opt)

            # take additional samples if required
            flag = update_samples(estimator, n_opt)
            flag && break

        else

            while varest(estimator) > θ*ϵ^2
                max_index = find_index_with_max_var_over_cost(estimator) 

                # increase the number of samples already taken
                n_opt = next_number_of_samples(estimator, max_index)

                # print optimal number of samples
                estimator[:verbose] && print_optimal_nb_of_samples(estimator, n_opt)

                # take additional samples
                flag = update_samples(estimator, n_opt)
                flag && break

                # recompute splitting parameter
                θ = T1 <: SL ? 0.5 : estimator[:do_mse_splitting] ? compute_splitting(estimator, ϵ) : estimator[:min_splitting]

                # check next iteration
                estimator[:verbose] && print_qmc_convergence(estimator, ϵ, θ)
            end

            flag && break

        end

        # check for convergence
        if T1 <: SL
            is_converged = true
        else
            # show status
            estimator[:verbose] && print_status(estimator)
            estimator[:verbose] && L ≥ 2 && print_mse_analysis(estimator, ϵ, θ)

            # check if the new level exceeds the maximum level
			is_converged = L > estimator[:min_index_set_param] && converged(estimator, ϵ, θ)
            if !is_converged && max_level_exceeded(estimator) 
                estimator[:verbose] && warn_max_level(estimator)
                is_converged = true
            else
                # update level
                L += 1
            end
        end
    end

    if !flag
        # update boundary in case of AD
        T1 <: AD && update_boundary(estimator)

        # print convergence status
        estimator[:verbose] && print_convergence(estimator, T1 <: SL ? true : converged(estimator, ϵ, θ))
    end
end

## Unbiased estimator routine ##
function _run(estimator::Estimator{<:U, <:AbstractSampleMethod}, ϵ::Real)

    # print status
    estimator[:verbose] && print_header(estimator, ϵ)

    # update index set
    index_set = collect(keys(pmf(estimator)))
    estimator[:verbose] && print_index_set(estimator, index_set)

    # add new indices to the index set
    for index in index_set
        index ∈ current_index_set(estimator) || push!(estimator, index)
    end

    # main loop
    is_converged = varest(estimator) ≤ ϵ^2
    flag = false
    while !is_converged

        # print status
        estimator[:verbose] && print_status(estimator)

        # increase the number of samples already taken
        n_opt = next_number_of_samples(estimator)

        # print optimal number of samples
        estimator[:verbose] && print_optimal_nb_of_samples(estimator, n_opt)

        # take additional samples
        flag = update_samples(estimator, n_opt)
        flag && break

        # update pmf
        update_pmf(estimator)
        estimator[:verbose] && print_pmf(estimator)

        # check next iteration
        estimator[:verbose] && print_unbiased_convergence(estimator, ϵ)

        is_converged = varest(estimator) ≤ ϵ^2

        # print index set
        estimator[:verbose] && print_index_set(estimator, index_set)
    end

    if !flag
        # print convergence status
        estimator[:verbose] && print_convergence(estimator, true)
    end
end
