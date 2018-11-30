## run.jl : entry point method for simulating an Estimator
#
# Runs the Estimator for the given tolerance, and returns a History object.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

"""
```julia
run(estimator, ϵ)
```
Run the estimator and compute the expected value of the quantity of interest up to the given tolerance(s) `ϵ`.

# Examples
```jldoctest
julia>

```
"""
function run(estimator::Estimator, tols::AbstractVector{<:Real})

    # input checking
    all(tols.>0) || throw(ArgumentError(string("supplied tolerance(s) must be positive, got ", tols)))

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

    # start level is 0
    L = 0

    # MSE splitting parameter
    θ = min_splitting(estimator)

    # main loop
    while L ≤ 2 || !converged(estimator, ϵ, θ)

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

        # value of the MSE splitting parameter
        θ = do_mse_splitting(estimator) ? compute_splitting(estimator, ϵ) : min_splitting(estimator)

        # evaluate optimal number of samples
        n_opt = Dict(τ => optimal_nb_of_samples(estimator, τ, ϵ, θ) for τ in keys(estimator))

        # print optimal number of samples
        verbose(estimator) && print_optimal_nb_of_samples(estimator, n_opt)

        # take additional samples
        for τ in keys(estimator)
            n_due = n_opt[τ] - nb_of_samples(estimator, τ)
            n_due > 0 && sample!(estimator, τ, n_due)
        end

        if T <: SL
            break
        else
            # show status
            verbose(estimator) && print_status(estimator)
            verbose(estimator) && L ≥ 2 && print_mse_analysis(estimator, ϵ, θ)

            # check if the new level exceeds the maximum level
            if !converged(estimator, ϵ, θ) && ( sz(estimator) ≥ max_index_set_param(estimator) ) 
                verbose(estimator) && warn_max_level(estimator)
                break
            end

            # update level
            L += 1
        end
    end

    # print convergence status
    verbose(estimator) && print_convergence(estimator, converged(estimator, ϵ, θ))
end

## Main routine ##
function _run(estimator::Estimator{<:MG, <:MC}, N::Integer)
    # print status
    verbose(estimator) && print_header(estimator, ϵ)

    # steps
    step = 1000

    # push all levels to the index set
    for index in get_index_set(estimator, mx_index_set_param(index_set))
        push!(estimator, index)
    end

    # steps of 1000
    for m in 1:div(N, step)+1
        n = N > step*wm ? step : N % step 

        # print status
        verbose(estimator) && print_status(estimator)

        # compute r param
        r = 1/2*(β(estimator) + γ(estimator))
        r[broadcast(|,isnan.(r),r.<=0)] = 1.5 # replace NaN's
        
        # compute indices where to take samples
        n_opt = Dict(i=>0 for i in keys(estimator))
        while sum(values(n_opt)) < n
            idx = Index(floor.(Int, randexpr.(log(2)*r))...)
            if idx ∉ keys(estimator)
                i_nearest = indmin([sum(abs.(idx_.-idx)) for idx_ in keys(estimator)])
                idx = index_set[i_nearest]
            end
            n_opt[idx] += 1
        end

        for τ in sort(keys(estimator))
            n_opt[τ] > 0 && sample!(estimator, τ, n_opt[τ])
        end

        # show status
        estimator.verbose && print_mse_analysis(estimator,ϵ,θ)
    end
end
