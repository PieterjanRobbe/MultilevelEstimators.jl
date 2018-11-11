## multiindex_quasi_monte_carlo.jl : run Multi-Index Quasi-Monte Carlo estimator

## main routine ##
function _run(estimator::MultiIndexQuasiMonteCarloEstimator, ϵ::T where {T<:Real})

    # print status
    estimator.verbose && print_header(estimator,ϵ)

    # index size parameter is 0
    level = 0

    # loop variables
    converged = false

    # main loop
    while !converged 

        # update index set
        index_set = new_index_set(estimator,level)

        # obtain initial variance estimate
        N0 = estimator.nb_of_warm_up_samples
        for index in index_set
            if !haskey(estimator.samples[1],index)
                sample!(estimator,index,N0)
            end
        end

        # add new level to the index set
        for index in index_set
            push!(estimator,index)
        end

        # print status
        estimator.verbose && print_status(estimator)

        # value of the MSE splitting parameter
        θ = estimator.do_splitting ? compute_splitting(estimator,ϵ) : 1/2

        while varest(estimator) > θ*ϵ^2 # doubling algorithm

            # find index with maximum variance
            idx = argmax(varest.(estimator, keys(estimator))./cost.(estimator,keys(estimator)))
            max_index = keys(estimator)[idx] 

            # increase the number of samples already taken
            if estimator.sample_multiplication_factor == 2 # (round to nearest power of two)
                n_opt = nextpow2(estimator.nsamples[max_index]+1) # + 1 to increase amount
            elseif estimator.sample_multiplication_factor <= 1
                n_opt = estimator.nsamples[max_index] + 1 # add 1 sample
            else
                n_opt = ceil(Int,estimator.nsamples[max_index]*estimator.sample_multiplication_factor)
            end

            # print optimal number of samples
            estimator.verbose && print_number_of_samples(estimator,Dict(max_index=>n_opt))

            # take additional samples
            n_due = n_opt - estimator.nsamples[max_index]
            n_due > 0 && sample!(estimator,max_index,n_due)

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
        level = level + 1

        # check if the new level exceeds the maximum level
        if max_level_exceeded(estimator,level,converged)
            estimator.verbose && warn_max_level(estimator)
            break
        end
    end

    # update maximum active set
    update_max_active(estimator)

    # print convergence status
    estimator.verbose && print_convergence(estimator,converged)
end
