## multilevel_monte_carlo.jl : run Monte Carlo estimator

## TODO :: THREE EXTENSIONS
# 1. splitting
# 2. do_regression (avoid expensive warm_up_samples, works hand in hand with 3.)
# 3. repeat sampling until variance is smaller than 1.1*θ*ϵ^2 (avoids adding extra level unnecessary)
# 4. samples0 ?
############## 5. max_index !!!
############## 6. in continuate, should only use set of keys that is currently in use!
##############    FIX: append current_index_set to Estimator
##############    let keys(estimator) return only those keys
############## 7. implement better bias formula
############## 8. PROBLEM ::: eigenfunctions are different !!!!!!!!
##############    check this with original implementation ???
##############    ===> SeparableCovarianceFunction!!!!!!!!
##############    NO! THIS IS NO SOLUTION BECAUSE OF SMOOTHNESS!!!
# 9. add possibility for own parallel_sample! function into estimator
## main routine ##
# 10. MLMC with multiple
############## 11. add rates
# 12. catch empty bias / variance at start of loop
# 13. merge with MonteCarloEstimator
function _run(estimator::MultiLevelMonteCarloEstimator, ϵ::T where {T<:Real})

    # print status
    estimator.verbose && print_header(estimator,ϵ)

    # start level is Level(0)
    level = Level(0)

    # loop variables
    converged = false

    # to avoid error about empty set 
    ########################
    while isempty(keys(estimator)) || ( level <= (2,) ) || !converged 

        # add new level to the index set
        push!(estimator,level)

        # value of the MSE splitting paramter
        θ = 0.5 # TODO catch empty bias/est here

        # obtain initial variance estimate
        N0 = estimator.nb_of_warm_up_samples
        if !haskey(estimator.samples,level) # if no samples have been taken yet on this level
            sample!(estimator,level,N0,true)
        end

        estimator.verbose && print_status(estimator)

        # evaluate optimal number of samples
        n_opt = Dict{Index,Int}()
        all_sum = sum(sqrt.([var(estimator,level)*cost(estimator,level) for level in keys(estimator)])) 
        for tau in keys(estimator)
            n_opt[tau] = ceil.(Int,1/(θ*ϵ^2) * sqrt(var(estimator,tau)/cost(estimator,tau)) * all_sum)
        end

        estimator.verbose && print_number_of_samples(estimator,n_opt)

        # take additional samples
        for tau in keys(estimator)
            n_due = n_opt[tau] - estimator.nsamples[tau]
            n_due > 0 && sample!(estimator,tau,n_due,false)
        end

        # show status
        estimator.verbose && print_mse_analysis(estimator,ϵ)

        # check convergence
        converged = bias(estimator)^2 <= (1-θ)*ϵ^2

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

## Multilevel Monte Carlo parallel sampling ##
function parallel_sample!(estimator::MultiLevelMonteCarloEstimator,index::Index,istart::N,iend::N) where {N<:Integer}

    # parallel sampling
    wp = CachingPool(workers())
    f(i) = estimator.sample_function(index,get_point(estimator.number_generator,i),estimator.user_data)
    t = @elapsed all_samples = pmap(wp,f,istart:iend)

    # extract samples
    samples = last.(all_samples)
    dsamples = first.(all_samples)

    # append samples
    append!(estimator.samples[index],dsamples)
    estimator.nsamples[index] += iend-istart+1
    estimator.total_work[index] += estimator.use_cost_model ? (iend-istart+1)*estimator.cost_model(index) : t
end

## inspector functions ##
cost(estimator::MultiLevelMonteCarloEstimator,index::Index) = estimator.total_work[index]/estimator.nsamples[index]

function point_with_max_var(estimator::MultiLevelMonteCarloEstimator)
    v = zeros(estimator.nb_of_qoi)
    for index in keys(estimator)
        for i = 1:estimator.nb_of_qoi
            v[i] += var([estimator.samples[index][j][i] for j in 1:estimator.nsamples[index]])
        end
    end
    return indmax(v)
end

for f in [:mean :var :skewness]
    ex1 = :(
            function $(f)(estimator::MultiLevelMonteCarloEstimator,index::Index)
                idx = point_with_max_var(estimator)
                $(f)([estimator.samples[index][j][idx] for j in 1:estimator.nsamples[index]])
            end
           )
    eval(ex1)
    ex2 = :(
            $(f)(estimator::MultiLevelMonteCarloEstimator) = sum([$(f)(estimator,index) for index in keys(estimator)])
           )
    eval(ex2)
end

function varest(estimator::MultiLevelMonteCarloEstimator,index::Index)
    n = estimator.nsamples[index]
    var(estimator,index)/n 
end

varest(estimator::MultiLevelMonteCarloEstimator) = sum([varest(estimator,index) for index in keys(estimator)])

function moment(estimator::MultiLevelMonteCarloEstimator,k::N where {N<:Integer},index::Index)
    idx = point_with_max_var(estimator)
    moment([estimator.samples[index][j][idx] for j in 1:estimator.nsamples[index]],k)
end

moment(estimator::MultiLevelMonteCarloEstimator,k::N where {N<:Integer}) = sum([moment(estimator,k,index) for index in keys(estimator)])

function bias(estimator::MultiLevelMonteCarloEstimator)
    L = maximum(keys(estimator)).+1
    θ = α(estimator,both=true,conservative=estimator.conservative_bias_estimate)
    2^(θ[1]+L[1]*θ[2])
end

mse(estimator::MultiLevelMonteCarloEstimator) = varest(estimator) + bias(estimator)^2

rmse(estimator::MultiLevelMonteCarloEstimator) = sqrt(mse(estimator))

## rates ##
function α(estimator::MultiLevelMonteCarloEstimator; both=false, conservative=true)
    max_idx = maximum(keys(estimator))
    if max_idx < (2,)
        return both ? [NaN, NaN] : NaN
    else
        x_start = conservative ? 1 : max_idx[1]-1
        x = x_start:max_idx[1]
        y = zeros(size(x))
        for i in x
            y[i-x_start+1] = mean(estimator,(i,))
        end
        θ = straight_line_fit(x,log2.(abs.(y)))
        return both ? θ : max(0.5,-θ[2])
    end
end

function β(estimator::MultiLevelMonteCarloEstimator; both=false)
    max_idx = maximum(keys(estimator))
    if max_idx < (2,)
        return both ? [NaN, NaN] : NaN
    else
        x = 1:max_idx[1]
        y = zeros(size(x))
        for i in x
            y[i] = var(estimator,(i,))
        end
        θ = straight_line_fit(x,log2.(y))
        return both ? θ : -θ[2]
    end
end

function γ(estimator::MultiLevelMonteCarloEstimator)
    max_idx = maximum(keys(estimator))
    if max_idx < (2,)
        return NaN
    else
        x = 1:max_idx[1]
        y = zeros(size(x))
        for i in x
            y[i] = cost(estimator,(i,))
        end
        return straight_line_fit(x,log2.(y))[2]
    end
end

function straight_line_fit(x::AbstractVector,y::AbstractVector)
    X = ones(length(x),2)
    X[:,2] = x
    X\y
end
