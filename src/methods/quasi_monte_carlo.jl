## quasi_monte_carlo.jl : run Quasi-Monte Carlo estimator

## main routine ##
function _run(estimator::QuasiMonteCarloEstimator, ϵ::T where {T<:Real})

    # print status
    estimator.verbose && print_header(estimator,ϵ)

    # QMC simulation is performed at Level(0)
    level = Level(0)
    push!(estimator,level)

    converged = false

    # obtain initial variance estimate
    N0 = estimator.nb_of_warm_up_samples
    if !haskey(estimator.samples[1],level)
        sample!(estimator,level,N0)
    end

    # print status
    estimator.verbose && print_status(estimator)

    while varest(estimator) > ϵ^2/2 # doubling algorithm

        # increase the number of samples already taken
        if estimator.sample_multiplication_factor == 2 # (round to nearest power of two)
            n = nextpow2(estimator.nsamples[level]+1) # + 1 to increase amount
        elseif estimator.sample_multiplication_factor <= 1
            n = estimator.nsamples[level] + 1 # add 1 sample
        else
            n = ceil(Int,estimator.nsamples[level]*estimator.sample_multiplication_factor)
        end

        # print optimal number of samples
        estimator.verbose && print_number_of_samples(estimator,Dict(level=>n))

        # take additional samples
        n_due = n - estimator.nsamples[level]
        n_due > 0 && sample!(estimator,level,n_due)
    end

    # print convergence status
    estimator.verbose && print_convergence(estimator,true)
end
