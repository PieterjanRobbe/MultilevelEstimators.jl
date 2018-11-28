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
