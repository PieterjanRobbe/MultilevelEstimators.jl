## multigrid_multilevel_monte_carlo.jl : run Multigrid Multilevel Monte Carlo estimator

## main routine ##
function _run(estimator::MultipleSemiCoarsenedMultiGridMultiIndexMonteCarloEstimator, ϵ::T where {T<:Real})

    # print status
    estimator.verbose && print_header(estimator,ϵ)

    # start level is maximum level already taken
	# TODO: tidy up ==>
	level = isempty(estimator.samples[1]) ? Level(0) : Level(maximum(getindex.(collect(keys(estimator.samples[1])),1)))
	for ell in 0:level[1]
		push!(estimator,Level(ell))
	end

    # loop variables
    converged = false

    # main loop
    while !converged

        # obtain initial variance estimate
		if !haskey(estimator.samples[1],level)
        	sample!(estimator,level,1)
        	push!(estimator,level)
		end

        # print status
        estimator.verbose && print_status(estimator)

        # value of the MSE splitting parameter
        θ = estimator.do_splitting ? compute_splitting(estimator,ϵ) : 1/2

        while varest(estimator) > θ*ϵ^2 # doubling algorithm

            # increase the number of samples already taken
            N = length(estimator.samples[1][Level(0)])
            if estimator.sample_multiplication_factor == 2 # (round to nearest power of two)
                n_opt = nextpow2(N+1) # + 1 to increase amount
            elseif estimator.sample_multiplication_factor <= 1
                n_opt = N + 1 # add 1 sample
            else
                n_opt = ceil(Int,N*estimator.sample_multiplication_factor)
            end

            # print optimal number of samples
			estimator.verbose && print_number_of_samples(estimator,Dict(Level(0)=>n_opt))

            # take additional samples
            r = 1/2*(β(estimator) + γ(estimator))
            r = isnan(r) ? 1.5 : r
			N_add = min.(floor.(Int,randexpr(log(2)*r,n_opt-N)),level[1])
			N_sum = Int64[sum(N_add.==ℓ) for ℓ = 0:level[1]]
			N_diff = append!(-diff(N_sum),N_sum[end]) # subtract samples on finer levels
			for tau in sort(collect(keys(estimator.samples[1])))
                n_due = N_sum[tau[1]+1]
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
        converged = ( level >= (2,) ) && ( bias(estimator)^2 <= (1-θ)*ϵ^2 || mse(estimator) <= ϵ^2 )

        # increase level
        level = level .+ 1

        # check if the new level exceeds the maximum level
        if !converged && ( level > (estimator.max_level,) ) 
            estimator.verbose && warn_max_level(estimator)
            break
        end
    end

    # print convergence status
    estimator.verbose && print_convergence(estimator,converged)
end
