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
        #TODO no clear for MG case
        estimator.index_set isa MG || clear(estimator) # prepare new run
    end

    return history
end

run(estimator::Estimator, tol::Real) = run(estimator, get_tols(estimator, tol))

## Main routine ##
function _run(estimator::Estimator{T, <:MC}, ϵ::Real) where T<:AbstractIndexSet
    # print status
    verbose(estimator) && print_header(estimator, ϵ)

    # warm-up
    index_set = get_index_set(TD(2), 2)
    for index in index_set
        if !contains_samples_at_index(estimator, index)
            set_weight(estimator, index, 1.0)
            sample!(estimator, index, nb_of_warm_up_samples(estimator))
            push!(estimator, index)
        end
    end

    # print status
    verbose(estimator) && print_status(estimator)

    while varest(estimator) > ϵ^2
        max_index_set = get_index_set(estimator, max_index_set_param(estimator))

        # compute weights
        set_weights(estimator, max_index_set)

        # print weights
        verbose(estimator) && print_rate_r(estimator)
        verbose(estimator) && print_weights(estimator)

        # compute amount of samples to take
        n_have = length(first(accumulator(estimator)))
        if sample_mul_factor(estimator) == 2 # (round to nearest power of two)
            n_to_take = nextpow2(n_have + 1) - n_have
        elseif sample_mul_factor(estimator) ≤ 1
            n_to_take = 1
        else
            n_to_take = max(1, ceil(Integer, sample_mul_factor(estimator)*n_have)) - n_have
        end

        # compute indices where to take samples
        n_opt = Dict(i=>0 for i in max_index_set)
        if length(keys(estimator)) == 1 # catch level 0
            n_opt[first(max_index_set)] = n_to_take
        else
            n_to_compute = n_to_take
            my_r = r(estimator)
            while sum(values(n_opt)) < n_to_take
                lvls = map(i->floor.(Int, randexpr(my_r[i], n_to_compute)), 1:length(my_r)) 
                idcs = map(i->Index(getindex.(lvls, i)...), 1:n_to_compute)
                for idx in max_index_set
                    n_opt[idx] += sum(broadcast(i->i==idx, idcs))
                end
                n_to_compute = n_to_take - sum(values(n_opt))
            end
        end
        verbose(estimator) && _print_index_set(keys(filter(p -> p.second > 0, n_opt)))

        # substract samples that will be taken with reuse
        for index in max_index_set
            R = CartesianIndices(UnitRange.(zero(eltype(T)), index.-1))
            for I in R
                n_opt[Index(I)] = max(0, n_opt[Index(I)] -  n_opt[index])
            end
        end

        # make sure that at least 2 samples are taken
        for index in max_index_set
            if n_opt[index] == 1
                n_opt[index] += 1
            end
        end

        # print optimal number of samples
        verbose(estimator) && print_optimal_nb_of_samples(estimator, filter(p -> p.second > 0, n_opt))

        # take the samples
        for τ in max_index_set
            if n_opt[τ] > 0
                # add indices to the estimator
                for i in CartesianIndices(UnitRange.(zero(τ), τ))
                    if !contains_samples_at_index(estimator, Index(i))
                        add_index(estimator, Index(i))
                        push!(estimator, Index(i))
                    end
                end
                # sample!
                sample!(estimator, τ, n_opt[τ])
            end
        end

        # show status
        verbose(estimator) && print_status(estimator)
        verbose(estimator) && print_var_est(estimator, ϵ, 1)
    end

    # show status
    verbose(estimator) && print_status(estimator)
    verbose(estimator) && print_mse_analysis(estimator, ϵ, 1)

    # print convergence status
    verbose(estimator) && print_convergence(estimator, converged(estimator, ϵ, 1))
end
