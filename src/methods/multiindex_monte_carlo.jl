## multiindex_monte_carlo.jl : Multi-Index Monte Carlo method
#
# Implementation of the Multi-Index Monte Carlo (MIMC) method.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## rates ##
for (f, g, sgn) in zip((:α, :β, :γ), (:mean, :var, :cost), (-1, -1, 1))
	eval(
		 quote
			 $f(estimator::Estimator{<:AbstractMI}) = $sgn.*getindex.($(Symbol("rates_", f))(estimator, 1:ndims(estimator)), 2)
			 $(Symbol("rates_", f))(estimator::Estimator{<:AbstractMI}, dir::Integer) = $(Symbol("rates_", f))(estimator, (maximum(getindex.(keys(estimator), dir)) + 1) * ntuple(i->i==dir, ndims(estimator)), dir)
			 function $(Symbol("rates_", f))(estimator::Estimator{<:AbstractMI}, idx::Index, dir::Integer)
				 m = idx[dir] - 1
				 if m < 2
					 return (NaN, NaN)
				 else
					 x = m:-1:1
					 y = map(xᵢ -> $g(estimator, idx .- xᵢ.*ntuple(i->i==dir, ndims(estimator))), 1:m)
					 idcs = .!isnan.(y)
					 θ = interp1(view(x, idcs), log2.(abs.(view(y, idcs))))
					 return tuple(θ...)
				 end
			 end
		 end)
end

## regression ##
function interpd(estimator::Estimator{I}, f::Function) where I<:AbstractMI
	idx_set = all_keys(estimator)
	A = [getindex(idx_set[i], j) for i in 2:length(idx_set), j in 1:ndims(estimator)]
	A = hcat(ones(eltype(A), size(A, 1)), A)
	y = map(index->log2(f(estimator, index)), view(idx_set, 2:length(idx_set)))
	A\y
end

for (f, sym) in zip([:var, :cost], [:β, :γ])
	eval(
		 quote
			 function $(Symbol("_regress_", f))(estimator::Estimator{<:MI, <:MC}, index::Index)
				 p = broadcast(dir->$(Symbol("rates_", sym))(estimator, index, dir), 1:ndims(estimator))
				 estimates = broadcast(dir->2^(p[dir][1]+index[dir]*p[dir][2]), 1:ndims(estimator))
				 estimate = mean(filter(!isnan, estimates))
				 if isnan(estimate)
					 p = interpd(estimator, $sym)
					 return 2 .^(p[1]+sum(p[2:end].*index))
				 else
				   	 return estimate
				 end
			 end
		 end)
end

regress_var(estimator::Estimator{<:MI, <:MC}, index::Index) = _regress_var(estimator::Estimator{<:MI, <:MC}, index::Index)
regress_cost(estimator::Estimator{<:MI, <:MC}, index::Index) = cost_model(estimator) isa EmtyFunction ? _regress_cost(estimator::Estimator{<:MI, <:MC}, index::Index) : cost_model(estimator, index)

# TODO also for ML
function regress_nb_of_samples(estimator::Estimator{<:AbstractIndexSet, <:MC}, index_set::Set{Index}, ϵ::Real, θ::Real)
	if do_regression(estimator) && L > 2
		vars = Dict(index=>regress_var(estimator, index) for index in index_set)
		costs = Dict(index=>regress_cost(estimator, index) for index in index_set)
		Σ_estimate = Σ(estimator)
		for index in keys(vars)
			Σ_estimate += sqrt(vars[index] * costs[index])
		end
		ns = Dict(index=>max(2, min(optimal_nb_of_samples(ϵ, θ, vars[index], costs[index], Σ_estimate), nb_of_warm_up_samples(estimator))) for index in index_set)
	else
		ict(index=>nb_of_warm_up_samples(estimator) index in index_set)
	end
end


























## bias ##
bias(estimator::Estimator{<:AbstractMI}) = bias(estimator, L(estimator))
function bias(estimator::Estimator{<:AbstractMI}, max_L)
    y = [abs(sum(mean.(estimator, boundary(estimator, cntr)))) for cntr in 1:max_L]
    start_L = robustify_bias_estimate(estimator) ? 1 : max_L-1
    p = interp1(start_L:max_L, log2.(y))
    2^(p[1]+(max_L+1)*p[2])
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
