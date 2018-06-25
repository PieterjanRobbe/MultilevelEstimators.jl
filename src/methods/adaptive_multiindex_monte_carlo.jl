## adaptive_multiindex_monte_carlo.jl : run Adaptive Multi-Index Monte Carlo estimator

## main routine ##
function _run(estimator::AdaptiveMultiIndexMonteCarloEstimator, ϵ::T where {T<:Real})

    # print status
    estimator.verbose && print_header(estimator,ϵ)

    # index size parameter is 0
    level = 0

    # dimension of index set
    d = ndims(estimator)

    # loop variables
    converged = false

    # main loop
    while !converged 

        # find index with largest "profit" and add to old set; enlarge active set
        # TODO is this the onyl difference? ===> move AD to MIMC/MIQMC
        index_set = typeof(estimator.current_index_set)()
        if level == 0
            max_index = ntuple(i->0,d)
            push!(index_set,max_index)
        else
            temp_active_set = collect(active_set(estimator))
            idx = indmax(profit.(estimator,temp_active_set))
            max_index = temp_active_set[idx]
            push!(estimator.old_index_set,max_index)
            print_largest_profit(max_index)
            for k in 1:d
             new_index = max_index .+ unit(k,d)
                if is_admissable(estimator.old_index_set,new_index)
                    if new_index ∈ get_index_set(estimator.max_search_space,estimator.max_level) 
                        push!(index_set,new_index)
                    else
                        warn_spill_index(max_index)
                        push!(estimator.spill_index_set,max_index)
                    end
                end
            end
        end

        # obtain initial variance estimate
        N0 = estimator.nb_of_warm_up_samples
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

        # value of the MSE splitting parameter
        θ = estimator.do_splitting ? compute_splitting(estimator,ϵ) : 1/2

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
        #if !converged && ( level > estimator.max_level ) 
        #    estimator.verbose && warn_max_level(estimator)
        #    break
        #end
        if isempty(active_set(estimator))
            estimator.verbose && warn_max_level(estimator)
            break
        end
    end

    # update maximum active set
    empty!(estimator.max_active_set)
    union!(estimator.max_active_set,active_set(estimator))
    union!(estimator.max_active_set,estimator.spill_index_set)

    # print convergence status
    estimator.verbose && print_convergence(estimator,converged)
end

# TODO for AMIQM, realize that this is not the best approach to compute gains...
function profit(estimator::AdaptiveMultiIndexTypeEstimator,index::Index)
    mean(estimator,index)/sqrt(var(estimator,index)*cost(estimator,index))
end

function bias(estimator::AdaptiveMultiIndexTypeEstimator; use_maximum=false::Bool)
    if use_maximum
        boundary = estimator.max_active_set
    else
        boundary = union(estimator.spill_index_set,active_set(estimator))
    end
    return length(boundary) < 2 ? NaN : sum(mean.(estimator,collect(boundary)))
end

#=

## inspector functions ##
function bias(estimator::MultiIndexTypeEstimator; use_maximum=false::Bool)
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
    if length(x) > 1
        start_level = estimator.conservative_bias_estimate ? 1 : length(x) - 1
        θ = straight_line_fit(x[start_level:end],log2.(abs.(y[start_level:end])))
        return 2^(θ[1]+(x[end]+1)*θ[2])
    else
        return NaN
    end
end

## rates ##
α(estimator::MultiIndexTypeEstimator) = α.(estimator,1:ndims(estimator))

function α(estimator::MultiIndexTypeEstimator, dir::Int64; both=false)
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

β(estimator::MultiIndexTypeEstimator) = β.(estimator,1:ndims(estimator))

function β(estimator::MultiIndexTypeEstimator, dir::Int64; both=false, start_idx=Index(zeros(Int64,ndims(estimator))...)::Index)
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
end

γ(estimator::MultiIndexTypeEstimator) = γ.(estimator,1:ndims(estimator))

function γ(estimator::MultiIndexTypeEstimator, dir::Int64; both=false, start_idx=Index(zeros(Int64,ndims(estimator))...)::Index)
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

# regress on all_sum: \sum \sqrt{V_\ell*C_\ell} (to avoid duplication)
function regress_all_sum(estimator::Estimator,index_set)
    all_sum = 0.
    for index in index_set
        if !haskey(estimator.samples[1],index) # no samples on this level yet; do regression on var/cost
            all_sum += sqrt(var_regress(estimator,index)*cost_regress(estimator,index))
        else # already samples on this level; use measured var/cost
            all_sum += sqrt(var(estimator,index)*cost(estimator,index))
        end
    end
    return all_sum
end

# regress variance and cost
for i in zip(["var","cost"],["β","γ"])
    ex = :(
           # regress index based on available information in estimator
           function $(Symbol(i[1],"_regress"))(estimator,index)
               value = Float64[] # collect estimated contribution in each direction
               for dir in 1:ndims(estimator)
                   if index[dir] > 2 # if enough levels are available in direction dir
                       θ = $(Symbol(i[2]))(estimator,dir,both=true,start_idx=index)
                       push!(value,2^(θ[1]+index[dir]*θ[2]))
                   end
               end
               if isempty(value)
                   # multi-dimensional regression
                   idx_set = setdiff(collect(keys(estimator.samples[1])),[tuple(zeros(ndims(estimator))...)])
                   m = length(idx_set); n = ndims(estimator)+1
                   X = ones(m,n)
                   X[1:m,2:n] = [getindex(idx_set[i],j) for i in 1:m, j in 1:n-1]
                   y = log2.([$(Symbol(i[1]))(estimator,index) for index in idx_set])
                   θ = X\y
                   return 2.^(θ[1]+sum(θ[2:end].*index))
               else
                   return mean(value)
               end
           end
           )
    eval(ex)
end

# regression
function regress(estimator::MultiIndexMonteCarloEstimator,index::Index,all_sum::T,ϵ::T,θ::T) where {T<:Real}
    var_est = var_regress(estimator,index)
    cost_est = cost_regress(estimator,index)
    n_opt = ceil.(Int,1/(θ*ϵ^2) * sqrt(var_est/cost_est) * all_sum)
    max(2,min(n_opt,estimator.nb_of_warm_up_samples))
end

=#
