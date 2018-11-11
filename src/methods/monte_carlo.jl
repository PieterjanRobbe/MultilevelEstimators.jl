## monte_carlo.jl : run Monte Carlo estimator

## main routine ##
function _run(estimator::Estimator{<:SL, <:_MC}, ϵ::T where {T<:Real})

    # print status
    estimator.verbose && print_header(estimator,ϵ)

    # MC simulation is performed at Level(0)
    level = Level(0)
    push!(estimator,level)
    
    # obtain initial variance estimate
    N0 = estimator.nb_of_warm_up_samples
    if !haskey(estimator.samples[1],level)
        sample!(estimator,level,N0)
    end

    # print status
    estimator.verbose && print_status(estimator)

    # evaluate optimal number of samples
    n = ceil.(Int,2/ϵ^2 * var(estimator,level))

    # print optimal number of samples
    estimator.verbose && print_number_of_samples(estimator,Dict(level=>n))

    # take additional samples
    n_due = n - estimator.nsamples[level]
    n_due > 0 && sample!(estimator,level,n_due)

    # print convergence status
    estimator.verbose && print_convergence(estimator,true)
end

# bias
bias(estimator::Estimator{<:SL}) = 0.

# rates
α(estimator::Estimator{<:SL}) = NaN
β(estimator::Estimator{<:SL}) = NaN
γ(estimator::Estimator{<:SL}) = NaN
