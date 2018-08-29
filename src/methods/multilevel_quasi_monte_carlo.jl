## multilevel_quasi_monte_carlo.jl : run Multilevel Quasi-Monte Carlo estimator

## main routine ##
function _run(estimator::MultiLevelQuasiMonteCarloEstimator, ϵ::T where {T<:Real})

    # print status
    estimator.verbose && print_header(estimator,ϵ)

    # start level is Level(0)
    level = Level(0)

    # loop variables
    converged = false

    # main loop
    while !converged 

        # obtain initial variance estimate
        N0 = estimator.nb_of_warm_up_samples
        if !haskey(estimator.samples[1],level)
            sample!(estimator,level,N0)
        end

        # add new level to the index set
        push!(estimator,level)

        # print status
        estimator.verbose && print_status(estimator)

        # value of the MSE splitting parameter
        θ = estimator.do_splitting ? compute_splitting(estimator,ϵ) : 1/2

        while varest(estimator) > θ*ϵ^2 # doubling algorithm

            # find level with maximum variance
            max_level = (indmax([varest(estimator,tau)/cost(estimator,tau) for tau in keys(estimator)])-1,) # convert to level

            # increase the number of samples already taken
            if estimator.sample_multiplication_factor == 2 # (round to nearest power of two)
                n_opt = nextpow2(estimator.nsamples[max_level]+1) # + 1 to increase amount
            elseif estimator.sample_multiplication_factor <= 1
                n_opt = estimator.nsamples[max_level] + 1 # add 1 sample
            else
                n_opt = ceil(Int,estimator.nsamples[max_level]*estimator.sample_multiplication_factor)
            end

            # print optimal number of samples
            estimator.verbose && print_number_of_samples(estimator,Dict(max_level=>n_opt))

            # take additional samples
            n_due = n_opt - estimator.nsamples[max_level]
            n_due > 0 && sample!(estimator,max_level,n_due)

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

## Multilevel Quasi-Monte Carlo parallel sampling ##
function parallel_sample!(estimator::QuasiMonteCarloTypeEstimator,index::Index,istart::N,iend::N) where {N<:Integer}

    # parallel sampling
    wp = CachingPool(workers())
    f(i) = estimator.sample_function(index,get_point(estimator.number_generators[index],i[2],i[1]),estimator.user_data)
    nshifts = estimator.nb_of_shifts
    nqoi = estimator.nb_of_qoi
    t = @elapsed all_samples = pmap(wp,f,Base.Iterators.product(1:nshifts,istart:iend)) # first changes fastest

    # extract samples
    samples = last.(all_samples)
    dsamples = first.(all_samples)

    # append samples
    for j in 1:nshifts
        for i in 1:nqoi
            append!(estimator.samples[i,j][index],getindex.(dsamples[j,:],i))
            append!(estimator.samples0[i,j][index],getindex.(samples[j,:],i))
        end
    end
    estimator.nsamples[index] += iend-istart+1
    estimator.total_work[index] += estimator.use_cost_model ? (iend-istart+1)*estimator.cost_model(index) : t
end

## inspector functions ##
point_with_max_var_est(estimator::QuasiMonteCarloTypeEstimator) = indmax([var(sum([mean.([estimator.samples[q,s][idx] for s in 1:estimator.nb_of_shifts]) for idx in keys(estimator)])) for q in 1:estimator.nb_of_qoi])

for f in [:mean :var :skewness]
    ex = :(
           function $(f)(estimator::QuasiMonteCarloTypeEstimator,index::Index)
               idx = point_with_max_var_est(estimator)
               mean(Float64[$(f)(estimator.samples[idx,s][index]) for s in 1:estimator.nb_of_shifts])
           end
          )
    eval(ex)
end

for f in [:mean :var] # for samples0
    ex = :(
           function $(Symbol(f,0))(estimator::QuasiMonteCarloTypeEstimator,index::Index)
               idx = point_with_max_var_est(estimator)
               mean(Float64[$(f)(estimator.samples0[idx,s][index]) for s in 1:estimator.nb_of_shifts])
           end
         )
    eval(ex)
end

for f in [:mean :var :skewness]
    ex = :(
           function $(f)(estimator::QuasiMonteCarloTypeEstimator)
               idx = point_with_max_var_est(estimator)
               mean([sum([$(f)(estimator.samples[idx,s][index]) for index in keys(estimator)]) for s in 1:estimator.nb_of_shifts])
           end
          )
    eval(ex)
end

function varest(estimator::QuasiMonteCarloTypeEstimator,index::Index)
    idx = point_with_max_var_est(estimator)
    var(Float64[mean(estimator.samples[idx,s][index]) for s in 1:estimator.nb_of_shifts])
end

function varest(estimator::QuasiMonteCarloTypeEstimator)
    idx = point_with_max_var_est(estimator)
    var(sum(Float64[mean(estimator.samples[idx,s][index]) for s in 1:estimator.nb_of_shifts] for index in keys(estimator)))
end

function moment(estimator::QuasiMonteCarloTypeEstimator,k::N where {N<:Integer},index::Index)
    idx = point_with_max_var_est(estimator)
    mean(Float64[moment(estimator.samples[idx,s][index],k) for s in 1:estimator.nb_of_shifts])
end

function moment(estimator::QuasiMonteCarloTypeEstimator,k::N where {N<:Integer})
    idx = point_with_max_var_est(estimator)
    mean(sum(Float64[moment(estimator.samples[idx,s][index],k) for s in 1:estimator.nb_of_shifts] for index in keys(estimator)))
end
