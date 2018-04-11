## monte_carlo.jl : run Monte Carlo estimator

## main routine ##
function _run(estimator::MonteCarloEstimator, ϵ::T where {T<:Real})

    # print status
    estimator.verbose && print_header(estimator,ϵ)

    # MC simulation is performed at Level(0)
    level = Level(0)
    
    # obtain initial variance estimate
    N0 = estimator.nb_of_warm_up_samples
    if !haskey(estimator,level)
        sample!(estimator,level,N0,true)
    end

    estimator.verbose && print_status(estimator)

    # evaluate optimal number of samples
    n = ceil.(Int,2/ϵ^2 * var(estimator,level))

    estimator.verbose && print_number_of_samples(estimator,Dict(level=>n))

    # take additional samples
    n_due = n - estimator.nsamples[level]
    n_due > 0 && sample!(estimator,level,n_due,false)

    # print convergence status
    estimator.verbose && print_convergence(estimator)
end

## Monte Carlo parallel sampling ##
function parallel_sample!(estimator::MonteCarloEstimator,index::Index,istart::N,iend::N) where {N<:Integer}

    # parallel sampling
    wp = CachingPool(workers())
    f(i) = estimator.sample_function(index,get_point(estimator.number_generator,i),estimator.user_data)
    t = @elapsed samples = pmap(wp,f,istart:iend)

    # append samples
    append!(estimator.samples[index],samples)
    estimator.nsamples[index] += iend-istart+1
    estimator.total_work[index] += estimator.use_cost_model ? (iend-istart+1)*estimator.cost_model(index) : t
end

## inspector functions ##
function point_with_max_var(estimator::MonteCarloEstimator,index::Index)
    v = zeros(estimator.nb_of_qoi)
    for i = 1:estimator.nb_of_qoi
        v[i] = var([estimator.samples[index][j][i] for j in 1:estimator.nsamples[index]])
    end
    return indmax(v)
end

for f in [:mean :var :skewness]
    ex1 = :(
        function $(f)(estimator::MonteCarloEstimator,index::Index)
            idx = point_with_max_var(estimator,index)
            $(f)([estimator.samples[index][j][idx] for j in 1:estimator.nsamples[index]])
        end
    )
    eval(ex1)
    ex2 = :(
            $(f)(estimator::MonteCarloEstimator) = sum([$(f)(estimator,index) for index in keys(estimator)])
    )
    eval(ex2)
end

function varest(estimator::MonteCarloEstimator,index::Index)
    n = estimator.nsamples[index]
    var(estimator,index)/n 
end

varest(estimator::MonteCarloEstimator) = sum([varest(estimator,index) for index in keys(estimator)])

function moment(estimator::MonteCarloEstimator,k::N where {N<:Integer},index::Index)
    idx = point_with_max_var(estimator,index)
    moment([estimator.samples[index][j][idx] for j in 1:estimator.nsamples[index]],k)
end

moment(estimator::MonteCarloEstimator,k::N where {N<:Integer}) = sum([moment(estimator,k,index) for index in keys(estimator)])

mse(estimator::MonteCarloEstimator) = varest(estimator)

rmse(estimator::MonteCarloEstimator) = sqrt(mse(estimator))
