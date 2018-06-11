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
        all_sum = ( level > 2 && estimator.do_regression ) ? regress_all_sum(estimator,index_set) : 0.
        for index in index_set
            if !haskey(estimator.samples[1],index)
                N0_ = ( level > 2 && estimator.do_regression ) ? regress(estimator,index,all_sum,ϵ,θ) : N0
                sample!(estimator,index,N0_)
            end
        end

        # add new level to the index set
        for index in index_set
            push!(estimator,index)
        end

        # print status
        estimator.verbose && print_status(estimator)

        # value of the MSE splitting paramter
        θ = ( estimator.do_splitting && maximum(maximum(keys(estimator))) > 1 ) ? compute_splitting(estimator,ϵ) : 1/2

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
    
    # print convergence status
    estimator.verbose && print_convergence(estimator,converged)
end

# NEW APPROACH
# DO LINEAR FIT OF BIAS!!!
# TODO: check if this is efficient...
function bias(estimator::IndexTypeEstimator; use_maximum=false::Bool)
    index_set = use_maximum ? keys(estimator.samples[1]) : keys(estimator)
    x = Int64[]
    y = Float64[]
    level = 1
    max_reached = false
    while !max_reached
        boundary = setdiff(get_index_set(estimator.method,level),get_index_set(estimator.method,level-1))
        if isempty(∩(index_set,boundary))
            max_reached = true
        else
            push!(x,level)
            push!(y,sum(mean.(estimator,boundary)))
            level += 1
        end
    end
    println("----- bias -----")
    @show x
    @show y
    if length(x) > 1
        start_level = estimator.conservative_bias_estimate ? 1 : length(x) - 1
        θ = straight_line_fit(x[start_level:end],log2.(abs.(y[start_level:end])))
        return 2^(θ[1]+(x[end]+1)*θ[2])
    else
        return NaN
    end
end

#=
    @show boundary = use_maximum ? isempty(estimator.max_boundary) ? estimator.current_boundary : estimator.max_boundary : estimator.current_boundary
    B = 0. # bias contribution
    for index in boundary 
        value = Float64[] # collect estimated contribution in each direction
        for dir in 1:ndims(estimator)
            println("===============")
            @show index
            @show dir
            println("===============")
            if index[dir] > 2 # if enough levels are available in direction dir
                θ = α(estimator,dir,both=true,conservative=estimator.conservative_bias_estimate,start_idx=index)
                @show θ
                push!(value,2^(θ[1]+index[dir]*θ[2]))
            end
        end
        @show index
        @show value
        B += isempty(value) ? 0. : mean(value) # take mean value over all estimates
    end
    return iszero(B) ? NaN : B
end
=#

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
        return both ? θ : -θ[2]
    end
end

# TODO β.(estimator,1:ndims(estimator))
β(estimator::IndexTypeEstimator) = tuple([β(estimator,i) for i in 1:ndims(estimator)]...)

function β(estimator::IndexTypeEstimator, dir::Int64; both=false, start_idx=Index(zeros(Int64,ndims(estimator))...)::Index)
    max_idx = all(start_idx.==0) ? maximum(getindex.(keys(estimator),dir)) : start_idx[dir]-1
    if max_idx < 2
        return both ? [NaN, NaN] : NaN
    else
        x = 1:max_idx
        start_vec = [start_idx...]
        start_vec[dir] = 0
        new_start_idx = Index(start_vec...)
        y = [var(estimator,new_start_idx.+i.*unit(dir,ndims(estimator))) for i in x]
        θ = straight_line_fit(x,log2.(y))
        return both ? θ : -θ[2]
    end
#=
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
    =#
end

γ(estimator::IndexTypeEstimator) = tuple([γ(estimator,i) for i in 1:ndims(estimator)]...)

function γ(estimator::IndexTypeEstimator, dir::Int64; both=false, start_idx=Index(zeros(Int64,ndims(estimator))...)::Index)
    max_idx = all(start_idx.==0) ? maximum(getindex.(keys(estimator),dir)) : start_idx[dir]-1
    if max_idx < 2
        return both ? [NaN, NaN] : NaN
    else
        x = 1:max_idx
        start_vec = [start_idx...]
        start_vec[dir] = 0
        new_start_idx = Index(start_vec...)
        y = [cost(estimator,new_start_idx.+i.*unit(dir,ndims(estimator))) for i in x]
        θ = straight_line_fit(x,log2.(y))
        return both ? θ : θ[2]
    end
end

# compute optimal value of MSE splitting parameter
#function compute_splitting(estimator::IndexTypeEstimator,ϵ::T where {T<:Float64})
#    bias_est = bias(estimator,use_maximum=true)
#    compute_splitting(bias_est,ϵ)
#end

# regress on all_sum: \sum \sqrt{V_\ell*C_\ell}
function regress_all_sum(estimator::Estimator,index_set)
    all_sum = 0.
    for index in index_set
        @show index
        if !haskey(estimator.samples[1],index) # no samples on this level yet; do regression on var/cost
          @show  all_sum += sqrt(var_regress(estimator,index)*cost_regress(estimator,index))
        else # already samples on this level; use measured var/cost
          @show  all_sum += sqrt(var(estimator,index)*cost(estimator,index))
        end
    end
    return all_sum
end

# regress variance on index based on available information in estimator
function var_regress(estimator,index)
    value = Float64[] # collect estimated contribution in each direction
    for dir in 1:ndims(estimator)
        if index[dir] > 2 # if enough levels are available in direction dir
            θ = β(estimator,dir,both=true,start_idx=index)
            push!(value,2^(θ[1]+index[dir]*θ[2]))
        end
    end
    if isempty(value)
        # multi-dimensional regression
        idx_set = setdiff(collect(keys(estimator.samples[1])),[tuple(zeros(ndims(estimator))...)])
        m = length(idx_set); n = ndims(estimator)+1
        X = ones(m,n)
        X[1:m,2:n] = [getindex(idx_set[i],j) for i in 1:m, j in 1:n-1]
        y = [var(estimator,index) for index in idx_set]
        θ = X\y
        return 2.^(θ[1]+sum(θ[2:end].*index))
    else
        return mean(value)
    end
end

# regress cost on index based on available information in estimator
function cost_regress(estimator,index)
    value = Float64[] # collect estimated contribution in each direction
    for dir in 1:ndims(estimator)
        if index[dir] > 2 # if enough levels are available in direction dir
            θ = γ(estimator,dir,both=true,start_idx=index)
            push!(value,2^(θ[1]+index[dir]*θ[2]))
        end
    end
    if isempty(value)
        # multi-dimensional regression
        idx_set = setdiff(collect(keys(estimator.samples[1])),[tuple(zeros(ndims(estimator))...)])
        m = length(idx_set); n = ndims(estimator)+1
        X = ones(m,n)
        X[1:m,2:n] = [getindex(idx_set[i],j) for i in 1:m, j in 1:n-1]
        y = [cost(estimator,index) for index in idx_set]
        θ = X\y
        return 2.^(θ[1]+sum(θ[2:end].*index)) # TODO
    else
        return mean(value)
    end
end

# regression
function regress(estimator::MultiIndexMonteCarloEstimator,index::Index,all_sum::T,ϵ::T,θ::T) where {T<:Real}
    @show var_est = var_regress(estimator,index)
    @show cost_est = cost_regress(estimator,index)
    @show all_sum
    @show θ
    @show ϵ
    n_opt = ceil.(Int,1/(θ*ϵ^2) * sqrt(var_est/cost_est) * all_sum)
    max(2,min(n_opt,estimator.nb_of_warm_up_samples))
end

# TODO merge the next two functions?
# save the current boundary
#=
function update_boundary(estimator::IndexTypeEstimator,level::N where {N<:Integer})
    empty!(estimator.current_boundary)
    boundary = setdiff(get_index_set(estimator.method,level+1),get_index_set(estimator.method,level))
    for idx in boundary
        push!(estimator.current_boundary,idx)
    end
end

# save the max boundary
function update_max_boundary(estimator::IndexTypeEstimator,level::N where {N<:Integer})
    empty!(estimator.max_boundary)
    boundary = setdiff(get_index_set(estimator.method,level+1),get_index_set(estimator.method,level))
    for idx in boundary
        push!(estimator.max_boundary,idx)
    end
end
=#

# regression of optimal number of samples at unknown index
#=
function regress(estimator::MultiIndexMonteCarloEstimator,level::Index,ϵ::T,θ::T) where {T<:Real}

    # total variance on already known levels

    # total variance on 

    v = Float[]
    d = ndims(estimator)
    for dir in 1:d
        if index[dir] > 2
            x = 1:index[dir]-1
            y = [var(estimator,i.*unit(dir,d)) for i in x]
            θ = straight_line_fit(x,log2.(y))
            push!(v,2^(θ[1]+index[dir]*θ[2]))
        end
    end
    return isempty(v) ? estimator.nb_of_warm_up_samples : mean(v)
end


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
    
    p = β(estimator,both=true)
    var_est = 2^(p[1]+level[1]*p[2])
    p = γ(estimator,both=true)
    cost_est = 2^(p[1]+level[1]*p[2])
    all_sum = sum(sqrt.([var(estimator,level)*cost(estimator,level) for level in setdiff(keys(estimator),level)])) 
    all_sum += sqrt(var_est*cost_est)
    n_opt = ceil.(Int,1/(θ*ϵ^2) * sqrt(var_est/cost_est) * all_sum)
    max(2,min(n_opt,estimator.nb_of_warm_up_samples))
end
=#


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
