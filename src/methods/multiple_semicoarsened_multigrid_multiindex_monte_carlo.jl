## multiple_semicoarsened_multigrid_multiindex_monte_carlo.jl : run MSG Multi-Index Monte Carlo estimator

## main routine ##
function _run(estimator::MultipleSemiCoarsenedMultiGridMultiIndexMonteCarloEstimator, ϵ::T where {T<:Real})

    # print status
    estimator.verbose && print_header(estimator,ϵ)

    # start level is maximum level already taken
    indices = collect(keys(estimator.samples[1]))

    level = 0
    while isempty(setdiff(get_index_set(estimator.method,level),indices))
        level += 1
    end
    level = max(0, level-1)

    for ℓ in indices
        push!(estimator,ℓ)
    end

    # loop variables
    converged = false

    # main loop
    while !converged

        # update index set
        # TODO make sure that this is not level + 1
        # is collect necessary???
        @show index_set = collect(new_index_set(estimator,level))


        # obtain initial variance estimate
        for index in index_set
            if !haskey(estimator.samples[1],index)
                sample!(estimator,index,1)
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
            @show r[isnan.(r)] = 1.5 # replace NaN's

            N_sum = fill(0,1.+maximum.(collect(getindex.(index_set,j) for j in 1:ndims(estimator)))...)
            while sum(N_sum) < n_opt-N
                idx = Index(floor.(Int,randexpr.(log(2)*r))...)
                if idx ∈ index_set
                    N_sum[idx.+1...] += 1
                else
                    # TODO: test dit ???
                    i_nearest = indmin(abs.([sum(i.-idx) for i in index_set]))
                    idx_nearest = index_set[i_nearest]
                    N_sum[idx_nearest.+1...] += 1
                end
            end
            @show N_sum
            # difference

           ####### TODO fix this and done :  @show N_diff



            # TODO: formula for N ?
            # TODO: assume product structure 
            # TODO: ∏ᵢ₌₁ᵈ 2ᵅ⁽ⁱ⁾ˡ⁽ⁱ⁾ ... uitschrijven ...
            #@show r = 1/2*(β(estimator) + γ(estimator))
            #r[isnan.(r)] = 1.5 # replace NaN's
            # TODO: maybe do this adaptively? add more level samples to N_add if they are in the index set
            # until we reach n_opt-N
            # TODO BUT these samples are indices??? how to generate??? 
            #N_add = min.(floor.(Int,randexpr(log(2)*r,n_opt-N)),level[1])
            #N_sum = Int64[sum(N_add.==ℓ) for τ in estimator.current_index_set]
            # TODO should be Multi-Index difference...
            #N_diff = append!(-diff(N_sum),N_sum[end]) # subtract samples on finer levels
            # END TODO:

            for tau in sort(index_set)
                n_due = N_sum[tau.+1...]
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
    end

    # update maximum active set
    update_max_active(estimator)

    # print convergence status
    estimator.verbose && print_convergence(estimator,converged)
end
