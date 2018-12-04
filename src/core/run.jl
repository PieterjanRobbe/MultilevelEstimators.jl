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
	L = T <: MG ? max_sz(estimator) : 0
	set_sz(estimator, L)

	# add indices in MG case
	for ℓ in 0:L-1
		index_set = boundary(estimator, ℓ)
		for index in index_set
			push!(estimator, index)
		end
	end


	# MSE splitting parameter
	θ = min_splitting(estimator)

	#
	# main loop
	#
	while L ≤ 2 || !converged(estimator, ϵ, θ)

		# print level
		verbose(estimator) && print_level(estimator, L)

		# update index set
		index_set = boundary(estimator, L)
		verbose(estimator) && print_index_set(estimator, index_set)

		# update weights in unbiased MG estimation
		if T <: MG
			set_weights(estimator, index_set)
		end

		#
		# take warm-up samples if requires to obtain initial variance estimate
		#
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
		
		#
		# update samples to reduce variance
		#
		if T <: MG
			while varest(estimator)*θ > ϵ^2
				# compute weights
				set_weights(estimator)

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
				n_opt = Dict(i=>0 for i in keys(estimator))
				while sum(values(n_opt)) < n_to_take
					idx = Index(floor.(Int, randexpr.(r(estimator)))...)
					if idx ∈ keys(estimator)
						n_opt[idx] += 1
					end
				end

				# substract samples that will be taken with reuse
				for index in keys(estimator)
					R = CartesianIndices(UnitRange.(zero(eltype(T)), index.-1))
					for I in R
						n_opt[Index(I)] = max(0, n_opt[Index(I)] -  n_opt[index])
					end
				end

				# print optimal number of samples
				verbose(estimator) && print_optimal_nb_of_samples(estimator, n_opt)

				# take the samples
				for τ in keys(estimator)
					n_opt[τ] > 0 && sample!(estimator, τ, n_opt[τ])
				end

				# update splitting parameter
				θ = do_mse_splitting(estimator) ? compute_splitting(estimator, ϵ) : min_splitting(estimator)
			
				# show status
				verbose(estimator) && print_status(estimator)
				verbose(estimator) && L ≥ 2 && print_mse_analysis(estimator, ϵ, θ)
			end
		else
			# evaluate optimal number of samples
			n_opt = Dict(τ => optimal_nb_of_samples(estimator, τ, ϵ, θ) for τ in keys(estimator))

			# print optimal number of samples
			verbose(estimator) && print_optimal_nb_of_samples(estimator, n_opt)

			# take additional samples
			for τ in keys(estimator)
				n_due = n_opt[τ] - nb_of_samples(estimator, τ)
				n_due > 0 && sample!(estimator, τ, n_due)
			end
		end

		#
		# decide if additional levels should be added
		#
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
