## multigrid_multilevel_monte_carlo.jl : run Multigrid Multilevel Monte Carlo estimator

## main routine ##
function _run(estimator::MultiGridMultiLevelMonteCarloEstimator, ϵ::T where {T<:Real})

    # print status
    estimator.verbose && print_header(estimator,ϵ)

    # start level is 0
    level = Level(0)

    # loop variables
    converged = false

    # main loop
    while !converged

        # obtain initial variance estimate
        #N0 = estimator.nb_of_warm_up_samples
        # TODO: To determine N, I need to compute the variance ==> regress expected value 
        # TODO: NO! I think you can avooid this, but you need to regress the total variance of the estimator to compute splitting/varest
        #if !haskey(estimator.samples[1],level)
        #    N0_ = ( level > (2,) && estimator.do_regression ) ? regress(estimator,level,ϵ,θ) : N0 # regression
        #    sample!(estimator,level,N0_)
        #end

        # add new level to the index set
        push!(estimator,level)

        # print status
        estimator.verbose && print_status(estimator)

        # value of the MSE splitting parameter
        θ = estimator.do_splitting ? compute_splitting(estimator,ϵ) : 1/2

        while varest(estimator) > θ*ϵ^2 # doubling algorithm

            # increase the number of samples already taken
            if estimator.sample_multiplication_factor == 2 # (round to nearest power of two)
                n_opt = nextpow2(estimator.nsamples[(0,)]+1) # + 1 to increase amount
            elseif estimator.sample_multiplication_factor <= 1
                n_opt = estimator.nsamples[(0,)] + 1 # add 1 sample
            else
                n_opt = ceil(Int,estimator.nsamples[(0,)]*estimator.sample_multiplication_factor)
            end

            # print optimal number of samples
            estimator.verbose && print_number_of_samples(estimator,n_opt)

            # take additional samples
            r = 1/2*(β(estimator) + γ(estimator))
            r = isnan(r) ? 1.5 : r
            N = min.(floor.(Int,randexpr(r,n_opt)),level[1])
            N_sum = Int64[sum(N.==ℓ) for ℓ = 0:maximum(N)]
            for tau in keys(estimator)
                n_due = N_sum[tau[1]] - estimator.nsamples[tau]
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

# sample from exponential distribution with rate r
randexpr(r::Number,kwargs...) = randexp(kwargs...)/r

## Multilevel Monte Carlo parallel sampling ##
function parallel_sample!(estimator::MonteCarloTypeEstimator,index::Index,istart::N,iend::N) where {N<:Integer}

    # parallel sampling
    wp = CachingPool(workers())
    f(i) = estimator.sample_function(index,get_point(estimator.number_generators[index],i),estimator.user_data)
    t = @elapsed all_samples = pmap(wp,f,istart:iend)

    # extract samples
    samples = last.(all_samples)
    dsamples = first.(all_samples)

    # append samples
    for n_qoi in 1:estimator.nb_of_qoi
        append!(estimator.samples[n_qoi][index],getindex.(dsamples,n_qoi))
        append!(estimator.samples0[n_qoi][index],getindex.(samples,n_qoi))
    end
    estimator.nsamples[index] += iend-istart+1
    estimator.total_work[index] += estimator.use_cost_model ? (iend-istart+1)*estimator.cost_model(index) : t
end

## inspector functions ##
cost(estimator::MonteCarloTypeEstimator,index::Index) = estimator.total_work[index]/estimator.nsamples[index]

point_with_max_var(estimator::MonteCarloTypeEstimator) = indmax([ sum([var(estimator.samples[i][idx]) for idx in keys(estimator)]) for i in 1:estimator.nb_of_qoi ])

for f in [:mean :var :skewness]
    ex = :(
           function $(f)(estimator::MonteCarloTypeEstimator,index::Index)
               idx = point_with_max_var(estimator)
               $(f)(estimator.samples[idx][index])
           end
          )
    eval(ex)
end

for f in [:mean :var] # for samples0
    ex = :(
           function $(Symbol(f,0))(estimator::MonteCarloTypeEstimator,index::Index)
               idx = point_with_max_var(estimator)
               $(f)(estimator.samples0[idx][index])
           end
          )
    eval(ex)
end

for f in [:mean :var :skewness]
    ex = :(
           function $(f)(estimator::MonteCarloTypeEstimator)
               idx = point_with_max_var(estimator)
               sum(Float64[$(f)(estimator.samples[idx][index]) for index in keys(estimator)])
           end
          )
    eval(ex)
end

function varest(estimator::MonteCarloTypeEstimator,index::Index)
    n = estimator.nsamples[index]
    var(estimator,index)/n 
end

function varest(estimator::MonteCarloTypeEstimator)
    idx = point_with_max_var(estimator)
    sum([var(estimator.samples[idx][index])/estimator.nsamples[index] for index in keys(estimator)])
end

function moment(estimator::MonteCarloTypeEstimator,k::N where {N<:Integer},index::Index)
    idx = point_with_max_var(estimator)
    moment(estimator.samples[idx][index],k)
end

function moment(estimator::MonteCarloTypeEstimator,k::N where {N<:Integer})
    idx = point_with_max_var(estimator)
    sum([moment(estimator.samples[idx][index],k) for index in keys(estimator)])
end

function bias(estimator::MultiLevelTypeEstimator; use_maximum=false::Bool)
    max_idx = use_maximum ? maximum(keys(estimator.samples[1])) : maximum(keys(estimator))
    θ = α(estimator,both=true,conservative=estimator.conservative_bias_estimate,max_idx=max_idx)
    2^(θ[1]+(max_idx[1]+1)*θ[2])
end

mse(estimator::Estimator) = varest(estimator) + bias(estimator)^2

rmse(estimator::Estimator) = sqrt(mse(estimator))

## rates ##
function α(estimator::MultiLevelTypeEstimator; both=false, conservative=true, max_idx=maximum(keys(estimator))::Level) # optional arguments account for regression, MSE splitting etc.
    if max_idx < (2,)
        return both ? [NaN, NaN] : NaN
    else
        x_start = conservative ? 1 : max_idx[1]-1
        x = x_start:max_idx[1]
        y = [mean(estimator,(i,)) for i in x]
        θ = straight_line_fit(x,log2.(abs.(y)))
        return both ? θ : -θ[2]
    end
end

function β(estimator::MultiLevelTypeEstimator; both=false)
    max_idx = maximum(keys(estimator))
    if max_idx < (2,)
        return both ? [NaN, NaN] : NaN
    else
        x = 1:max_idx[1]
        y = Float64[]
        for i in x
            v = var(estimator,(i,))
            if !isnan(v)
                push!(y,v)
            end
        end
        x = x[1:length(y)]
        θ = straight_line_fit(x,log2.(y))
        if length(x) < 2
            return both ? [NaN, NaN] : NaN
        else
            return both ? θ : -θ[2]
        end
    end
end

function γ(estimator::MultiLevelTypeEstimator; both=false)
    max_idx = maximum(keys(estimator))
    if max_idx < (2,)
        return both ? [NaN, NaN] : NaN
    else
        x = 1:max_idx[1]
        y = [cost(estimator,(i,)) for i in x]
        θ = straight_line_fit(x,log2.(y))
        return both ? θ : θ[2]
    end
end

function straight_line_fit(x::AbstractVector,y::AbstractVector)
    X = ones(length(x),2)
    X[:,2] = x
    X\y
end

# regression of optimal number of samples at unknown level
function regress(estimator::MultiLevelMonteCarloEstimator,level::Level,ϵ::T,θ::T) where {T<:Real}
    p = β(estimator,both=true)
    var_est = 2^(p[1]+level[1]*p[2])
    p = γ(estimator,both=true)
    cost_est = 2^(p[1]+level[1]*p[2])
    all_sum = sum(sqrt.([var(estimator,level)*cost(estimator,level) for level in setdiff(keys(estimator),level)])) 
    all_sum += sqrt(var_est*cost_est)
    n_opt = ceil.(Int,1/(θ*ϵ^2) * sqrt(var_est/cost_est) * all_sum)
    max(2,min(n_opt,estimator.nb_of_warm_up_samples))
end

# compute optimal value of MSE splitting parameter
function compute_splitting(estimator::Estimator,ϵ::T where {T<:Real})
    bias_est = bias(estimator,use_maximum=true)
    isnan(bias_est) ? 0.5 : min(0.99, max(1/2,1-bias_est^2/ϵ^2))
end

