## multiindex_monte_carlo.jl : run Multi-Index Monte Carlo estimator

## main routine ##
function _run(estimator::MultiIndexMonteCarloEstimator, ϵ::T where {T<:Real})

    # print status
    estimator.verbose && print_header(estimator,ϵ)

    # index size parameter 0
    level = 0
    
    # loop variables
    converged = false

    # main loop
    while !converged 

        # obtain initial variance estimate
        N0 = estimator.nb_of_warm_up_samples
        index_set = get_index_set(estimator.method,level)
        for index in index_set
            if !haskey(estimator.samples[1],index)
                N0_ = ( level > 2 && estimator.do_regression ) ? regress(estimator,index,ϵ,θ) : N0 # regression
                sample!(estimator,index,N0_)
            end
        end

        # add new level to the index set
        for index in index_set
            push!(estimator,index)
        end
        update_boundary(estimator,level)

        # print status
        estimator.verbose && print_status(estimator)

        # value of the MSE splitting paramter
        θ = ( estimator.do_splitting && estimator.max_level_parameter.level > 1 ) ? compute_splitting(estimator,ϵ) : 1/2

        # evaluate optimal number of samples
        n_opt = Dict{Index,Int}()
        all_sum = sum(sqrt.([var(estimator,level)*cost(estimator,level) for level in keys(estimator)])) 
        for tau in keys(estimator)
            n_opt[tau] = ceil.(Int,1/(θ*ϵ^2) * sqrt(var(estimator,tau)/cost(estimator,tau)) * all_sum)
        end

        # print optimal number of samples
        estimator.verbose && print_number_of_samples(estimator,n_opt)

        # take additional samples
        for tau in keys(estimator)
            n_due = n_opt[tau] - estimator.nsamples[tau]
            n_due > 0 && sample!(estimator,tau,n_due)
        end

        # show status
        estimator.verbose && print_mse_analysis(estimator,ϵ,θ)

        # check convergence
        converged = ( level >= 2 ) && ( bias(estimator)^2 <= (1-θ)*ϵ^2 || mse(estimator) <= ϵ^2 )

        # increase level
        level = level + 1

        # check if the new level exceeds the maximum level
        if !converged && ( level > estimator.max_level ) 
            estimator.verbose && warn_max_level(estimator)
            break
        end
    end
    
    # save the current boundary for the next iterations
    update_max_boundary(estimator,level)

    # print convergence status
    estimator.verbose && print_convergence(estimator,converged)
end

function bias(estimator::IndexTypeEstimator; use_maximum=false::Bool)
    boundary = use_maximum ? estimator.max_boundary : estimator.current_boundary
    B = 0. # bias contribution
    for index in boundary 
        value = Float64[] # collect estimated contribution in each direction
        for dir in 1:ndims(estimator)
            if index[dir] > 2 # if enough levels are available in direction dir
                # fit from direction dir
                x_start = estimator.conservative_bias_estimate ? index[dir] - 1 : 2
                x = x_start:-1:1
                y = [mean(estimator,index.-i.*unit(dir,ndims(estimator))) for i in x]
                θ = straight_line_fit(x,log2.(abs.(y)))
                push!(value,2^(θ[1]))
            end
        end
        B += isempty(value) ? 0. : mean(value) # take mean value over all estimates
    end
    return iszero(B) ? NaN : B
end

## rates ##
α(estimator::IndexTypeEstimator) = tuple([α(estimator,i) for i in 1:ndims(estimator)]...)

function α(estimator::IndexTypeEstimator, dir::Int64; both=false)
    max_idx = maximum(getindex.(keys(estimator),dir))
    if max_idx < 2
        return both ? [NaN, NaN] : NaN
    else
        x = 1:max_idx
        y = [mean(estimator,i.*unit(dir,ndims(estimator))) for i in x]
        θ = straight_line_fit(x,log2.(abs.(y)))
        return both ? θ : max(0.5,-θ[2])
    end
end

β(estimator::IndexTypeEstimator) = tuple([β(estimator,i) for i in 1:ndims(estimator)]...)

function β(estimator::IndexTypeEstimator, dir::Int64; both=false)
    max_idx = maximum(getindex.(keys(estimator),dir))
    if max_idx < 2
        return both ? [NaN, NaN] : NaN
    else
        x = 1:max_idx
        y = Float64[]
        for i in x
            v = var(estimator,i.*unit(dir,ndims(estimator)))
            if !isnan(v)
                push!(y,v)
            end
        end
        x = x[1:length(y)]
        θ = straight_line_fit(x,log2.(y))
        return both ? θ : -θ[2]
    end
end

γ(estimator::IndexTypeEstimator) = tuple([γ(estimator,i) for i in 1:ndims(estimator)]...)

function γ(estimator::IndexTypeEstimator, dir::Int64; both=false)
    max_idx = maximum(getindex.(keys(estimator),dir))
    if max_idx < 2
        return both ? [NaN, NaN] : NaN
    else
        x = 1:max_idx
        y = [cost(estimator,i.*unit(dir,ndims(estimator))) for i in x]
        θ = straight_line_fit(x,log2.(y))
        return both ? θ : θ[2]
    end
end

# compute optimal value of MSE splitting parameter
function compute_splitting(estimator::IndexTypeEstimator,ϵ::T where {T<:Float64})
    bias_est = bias(estimator,use_maximum=true)
    compute_splitting(bias_est,ϵ)
end

# save the current boundary
function update_boundary(estimator::IndexTypeEstimator,level::N where {N<:Integer})
    empty!(estimator.current_boundary)
    boundary = setdiff(get_index_set(estimator.method,level),get_index_set(estimator.method,level+1))
    for idx in boundary
        push!(estimator.current_boundary,idx)
    end
end

# save the max boundary
function update_max_boundary(estimator::IndexTypeEstimator,level::N where {N<:Integer})
    empty!(estimator.max__boundary)
    boundary = setdiff(get_index_set(estimator.method,level),get_index_set(estimator.method,level+1))
    for idx in boundary
        push!(estimator.max_boundary,idx)
    end
end

# regression of optimal number of samples at unknown index
function regress(estimator::MultiIndexMonteCarloEstimator,level::Index,ϵ::T,θ::T) where {T<:Real}
    p = β(estimator,both=true)
    var_est = 2^(p[1]+level[1]*p[2])
    p = γ(estimator,both=true)
    cost_est = 2^(p[1]+level[1]*p[2])
    all_sum = sum(sqrt.([var(estimator,level)*cost(estimator,level) for level in setdiff(keys(estimator),level)])) 
    all_sum += sqrt(var_est*cost_est)
    n_opt = ceil.(Int,1/(θ*ϵ^2) * sqrt(var_est/cost_est) * all_sum)
    max(2,min(n_opt,estimator.nb_of_warm_up_samples))
end

#=

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

function bias(estimator::LevelTypeEstimator; max_idx=maximum(keys(estimator))::Level)
    L = max_idx.+1
    θ = α(estimator,both=true,conservative=estimator.conservative_bias_estimate,max_idx=max_idx)
    2^(θ[1]+L[1]*θ[2])
end

mse(estimator::Estimator) = varest(estimator) + bias(estimator)^2

rmse(estimator::Estimator) = sqrt(mse(estimator))

## rates ##
function α(estimator::LevelTypeEstimator; both=false, conservative=true, max_idx=maximum(keys(estimator))::Level) # optional arguments account for regression, MSE splitting etc.
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

function β(estimator::LevelTypeEstimator; both=false)
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

function γ(estimator::LevelTypeEstimator; both=false)
    max_idx = maximum(keys(estimator))
    if max_idx < (2,)
        return both ? [NaN, NaN] : NaN
    else
        x = 1:max_idx[1]
        y = zeros(size(x))
        for i in x
            y[i] = cost(estimator,(i,))
        end
        θ = straight_line_fit(x,log2.(y))
        return both ? θ : θ[2]
    end
end

function straight_line_fit(x::AbstractVector,y::AbstractVector)
    X = ones(length(x),2)
    X[:,2] = x
    X\y
end

# regression of optimal number of samples at unknown index
function regress(estimator::MultiLevelMonteCarloEstimator,level::Index,ϵ::T,θ::T) where {T<:Real}
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
function compute_splitting(estimator::LevelTypeEstimator,ϵ::T where {T<:Float64})
    L = max(maximum(keys(estimator.samples[1])),maximum(keys(estimator)))
    bias_est = bias(estimator,max_idx=L)
    min(0.99, max(1/2,1-bias_est^2/ϵ^2))
end
=#
