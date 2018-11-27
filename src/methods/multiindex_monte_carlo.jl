## multiindex_monte_carlo.jl : Multi-Index Monte Carlo method
#
# Implementation of the Multi-Index Monte Carlo (MIMC) method.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## main routine ##
function _run(estimator::Estimator{<:MI, <:MC}, ϵ::T where {T<:Real})

    # print status
    verbose(estimator) && print_header(estimator, ϵ)

    # index size parameter is 0
    set_L(estimator, 0)

    # loop variables
    θ = splitting(estimator)

    # main loop
    while L(estimator) ≤ 2 || !converged(estimator, ϵ, θ) 

        # update index set
        index_set = new_index_set(estimator, level)

        # obtain initial variance estimate
        for index in index_set
            if !contains_samples_at_index(estimator, index)
                if do_regression(estimator) && L(estimator) > 2
                    n = regress_nb_of_samples(estimator, index, ϵ, θ) # TODO !!! Σ contains unintialized indices
                    # must de regression on all indices! ==> regress(::Set{Index}) ?
                else
                    n = nb_of_warm_up_samples(estimator)
                end
                sample!(estimator, index, n)
            end
        end

        # add new level to the index set
        for index in index_set
            push!(estimator, index)
        end

        # print status
        verbose(estimator) && print_status(estimator)

        # value of the MSE splitting parameter
        θ = do_mse_splitting(estimator) ? compute_splitting(estimator, ϵ) : splitting(estimator)

        # evaluate optimal number of samples
        n_opt = Dict(τ => optimal_nb_of_samples(estimator, τ, ϵ, θ) for τ in keys(estimator))

        # print optimal number of samples
        verbose(estimator) && print_optimal_number_of_samples(estimator, n_opt)

        # take additional samples
        for τ in keys(estimator)
            n_due = n_opt[τ] - nb_of_samples(estimator, τ)
            n_due > 0 && sample!(estimator, τ, n_due)
        end

        # show status
        verbose(estimator) && print_status(estimator)
        verbose(estimator) && level ≥ 2 && print_mse_analysis(estimator, ϵ, θ)

        # increase level
        add_L(estimator)

        # check if the new level exceeds the maximum level
        if !converged(estimator, ϵ, θ) && max_level_exceeded(estimator, L(estimator))
            verbose(estimator) && warn_max_level(estimator)
            break
        end
    end

    # update max_L
    max_L(estimator) < L(estimator) && set_max_L(estimator, L(estimator)) 

    # update maximum active set
    update_max_active(estimator)

    # print convergence status
    verbose(estimator) && print_convergence(estimator, converged(estimator, ϵ, θ))
end

## bias ##
bias(estimator::Estimator{<:AbstractMI}) = bias(estimator, L(estimator))
function bias(estimator::Estimator{<:AbstractMI}, max_L)
    y = [abs(sum(mean.(estimator, boundary(estimator, cntr)))) for cntr in 1:max_L]
    start_L = robustify_bias_estimate(estimator) ? 1 : max_L-1
    p = interp1(start_L:max_L, log2.(y))
    2^(p[1]+(max_L+1)*p[2])
end
        

end

## inspector functions ##
function bias(estimator::MultiIndexTypeEstimator; use_maximum=false::Bool)


    if length(x) > 1
		if !estimator.conservative_bias_estimate
			return y[end]
		else
        start_level = estimator.conservative_bias_estimate ? 1 : length(x) - 1
        θ = straight_line_fit(x[start_level:end],log2.(abs.(y[start_level:end])))
        return 2^(θ[1]+(x[end]+1)*θ[2])
		end
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

## regression functions ##
function regress(estimator::MultiIndexTypeEstimator,index::Index,all_sum::T,ϵ::T,θ::T) where {T<:Real}
    var_est = var_regress(estimator,index)
    cost_est = cost_regress(estimator,index)
    n_opt = ceil.(Int,1/(θ*ϵ^2) * sqrt(var_est/cost_est) * all_sum)
    max(2,min(n_opt,estimator.nb_of_warm_up_samples))
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
                   return 2 .^(θ[1]+sum(θ[2:end].*index))
               else
                   return mean(value)
               end
           end
           )
    eval(ex)
end

## adaptivity ##
new_index_set(estimator::MultiIndexTypeEstimator, level::N where {N<:Integer}) = get_index_set(estimator.method,level)

function new_index_set(estimator::AdaptiveMultiIndexTypeEstimator, level::N where {N<:Integer})
    d = ndims(estimator) # dimension of index set
    # empty index set
    index_set = typeof(estimator.current_index_set)()
    if level == 0 # if first run
        max_index = ntuple(i->0,d)
        push!(index_set,max_index)
    else
        # find index with largest "profit" and add to old set; enlarge active set
        temp_active_set = collect(active_set(estimator))
        idx = argmax(profit.(estimator,temp_active_set))
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
    log_adaptive_index_set(estimator,keys(estimator) ∪ index_set,max_index)
    return index_set
end

# log indices (for plotting adaptive index set)
function log_adaptive_index_set(estimator,indices,max_index)
	d = ndims(estimator)
	dict = Dict{Index{d},Int64}()
	for index in indices
		if !haskey(dict,index)
			if index == max_index # max index
				dict[index] = 3
			elseif index ∈ estimator.spill_index_set # spill index
				dict[index] = 1
			elseif index ∈ estimator.old_index_set # old index
				dict[index] = 0
			else # active index
				dict[index] = 2
			end
		end
	end
	push!(estimator.adaptive_index_set,dict)
end

max_level_exceeded(estimator::MultiIndexTypeEstimator, level::N where {N<:Integer}, converged::Bool) = !converged && (level > estimator.max_level) 
max_level_exceeded(estimator::AdaptiveMultiIndexTypeEstimator, level::N where {N<:Integer}, converged::Bool) = !converged && isempty(active_set(estimator))

update_max_active(estimator::MultiIndexTypeEstimator) = nothing
function update_max_active(estimator::AdaptiveMultiIndexTypeEstimator)
    empty!(estimator.max_active_set)
    union!(estimator.max_active_set,active_set(estimator))
    union!(estimator.max_active_set,estimator.spill_index_set)
end

# TODO for AMIQM, realize that this is not the best approach to compute gains...
function profit(estimator::AdaptiveMultiIndexTypeEstimator,index::Index)
	if estimator.method isa MG
		abs(mean(estimator,index))
	else
		abs(mean(estimator,index))/sqrt(var(estimator,index)*cost(estimator,index))
	end
end

function bias(estimator::AdaptiveMultiIndexTypeEstimator; use_maximum=false::Bool)
    if use_maximum
        boundary = estimator.max_active_set
    else
        boundary = union(estimator.spill_index_set,active_set(estimator))
    end
	return length(boundary) < 2 ? NaN : abs.(sum(mean.(estimator,collect(boundary))))
end
