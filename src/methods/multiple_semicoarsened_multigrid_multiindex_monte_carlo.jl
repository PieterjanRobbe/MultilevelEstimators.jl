## multiple_semicoarsened_multigrid_multiindex_monte_carlo.jl : run MSG Multi-Index Monte Carlo estimator

## main routine ##
function _run(estimator::MultipleSemiCoarsenedMultiGridMultiIndexMonteCarloEstimator, ϵ::T where {T<:Real})

    # print status
    estimator.verbose && print_header(estimator,ϵ)

	# determine level parameter
	if estimator isa AdaptiveMultiIndexTypeEstimator
		level = length(estimator.old_index_set)
	else
		level = 0
		while isempty(setdiff(get_index_set(estimator.method,level),collect(keys(estimator.samples[1]))))
			level += 1
		end
		level = max(0, level-1)
	end
	index_set = level == 0 ? collect(new_index_set(estimator, 0)) : collect(keys(estimator))

    # loop variables
    converged = false

    # main loop
    while !converged

        # obtain initial variance estimate
        for index in index_set
            if !haskey(estimator.samples[1],index)
                sample!(estimator,index,2) ##### CHANGED 1 >> 2 ON 14/09 
                push!(estimator,index)
            end
        end

        # print status
        estimator.verbose && print_status(estimator)

        # value of the MSE splitting parameter
        θ = estimator.do_splitting ? compute_splitting(estimator,ϵ) : 1/2

        while varest(estimator) > θ*ϵ^2 # doubling algorithm

            # increase the number of samples already taken
            N = length(estimator.samples[1][zero_idx(estimator)])
            if estimator.sample_multiplication_factor == 2 # (round to nearest power of two)
                n_opt = nextpow2(N+1) # + 1 to increase amount
            elseif estimator.sample_multiplication_factor <= 1
                n_opt = N + 1 # add 1 sample
            else
                n_opt = ceil(Int,N*estimator.sample_multiplication_factor)
            end

            # print optimal number of samples
            estimator.verbose && print_number_of_samples(estimator,Dict(zero_idx(estimator)=>n_opt))

            # take additional samples
            r = 1/2*(β(estimator) + γ(estimator))
			r[broadcast(|,isnan.(r),r.<=0)] = 1.5 # replace NaN's
            N_sum = fill(0,1.+maximum.(collect(getindex.(index_set,j) for j in 1:ndims(estimator)))...)
			N_diff = copy(N_sum)
			while sum(N_diff) < n_opt-N
                idx = Index(floor.(Int,randexpr.(log(2)*r))...)
				if !(idx ∈ index_set)
					i_nearest = indmin([sum(abs.(i.-idx)) for i in index_set])
                    idx = index_set[i_nearest]
                end
                N_sum[idx.+1...] += 1
				N_diff[idx.+1...] += 1
				for i in Iterators.product(range.(1,idx)...)
					if i != idx && N_diff[i...] != 0
						N_diff[i...] -= 1
					end
				end
            end

            for tau in sort(index_set)
                n_due = N_diff[tau.+1...]
                n_due > 0 && sample!(estimator,tau,n_due)
            end

            # recompute splitting parameter
            θ = estimator.do_splitting ? compute_splitting(estimator,ϵ) : 1/2

            # check next iteration
            estimator.verbose && print_qmc_convergence(estimator,ϵ,θ)
        end

        # show status
        estimator.verbose && print_mse_analysis(estimator,ϵ,θ)

        # check convergence
        converged = ( level >= 2 ) && ( bias(estimator)^2 <= (1-θ)*ϵ^2 || mse(estimator) <= ϵ^2 )

        # increase level
        level += 1

        # check if the new level exceeds the maximum level
        if max_level_exceeded(estimator,level,converged)
            estimator.verbose && warn_max_level(estimator)
            break
        end

        # update index set
        index_set = converged ? index_set : collect(new_index_set(estimator,level))

    end

    # update maximum active set
    update_max_active(estimator)

    # print convergence status
    estimator.verbose && print_convergence(estimator,converged)

end
