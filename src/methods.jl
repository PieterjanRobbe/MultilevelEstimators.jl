## methods.jl : refactored methods from run.jl
#
# Implements refactored methods from run.jl.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2019

#
# basics
#
ndims(::Estimator{<:AbstractIndexSet{d}}) where d = d 

get_index_set(estimator::Estimator, sz) = get_index_set(estimator.index_set, sz)

get_tols(estimator::Estimator, tol::T) where T<:Real = estimator[:continuate] ? estimator[:continuation_mul_factor].^(estimator[:nb_of_tols]-1:-1:0)*tol : T[tol] 

mse(estimator::Estimator) = varest(estimator) + bias(estimator)^2

rmse(estimator::Estimator) = sqrt(mse(estimator))

converged(estimator::Estimator, ϵ::Real, θ::Real) = ( bias(estimator)^2 ≤ (1-θ)*ϵ^2 || mse(estimator) ≤ ϵ^2 )

max_level_exceeded(estimator::Estimator) = sz(estimator) ≥ estimator[:max_index_set_param]

#
# inspector functions: mean, var, varest...
#
qoi_with_max_var(estimator::Estimator{<:AbstractIndexSet, <:MC}) = argmax(map(n_qoi -> sum(var(samples_diff(estimator, n_qoi, index)) for index in keys(estimator)), 1:estimator[:nb_of_qoi]))

cost(estimator::Estimator, index::Index) = estimator[:cost_model] isa EmptyFunction ? time(estimator, index) : work(estimator, index)

for f in [:mean, :var]
    @eval begin

        $f(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index) = $f(samples_diff(estimator, qoi_with_max_var(estimator), index))

        $(Symbol(f, 0))(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index) = $f(samples(estimator, qoi_with_max_var(estimator), index))

        $f(estimator::Estimator{<:AbstractIndexSet, <:QMC}, index::Index) = mean($f(samples_diff(estimator, qoi_with_max_var(estimator), shift, index)) for shift in 1:estimator[:nb_of_shifts](index))

        $(Symbol(f, 0))(estimator::Estimator{<:AbstractIndexSet, <:QMC}, index::Index) = mean($f(samples(estimator, qoi_with_max_var(estimator), shift, index)) for shift in 1:estimator[:nb_of_shifts](index))

        $f(estimator::Estimator) = sum($f(estimator, index) for index in keys(estimator))

    end
end

varest(estimator::Estimator, index::Index) = var(estimator, index)/nb_of_samples(estimator, index)

varest(estimator::Estimator) =  sum(varest(estimator, index) for index in keys(estimator))

#
# rates 
#
for (f, g, sgn) in zip([:α, :β, :γ], [:mean, :var, :cost], [-1, -1, 1])
    @eval begin 

        $f(estimator::Estimator{<:SL}) = nothing

        $f(estimator::Estimator{<:AbstractIndexSet}) = $sgn.*getindex.(broadcast(i->$(Symbol("rates_", f))(estimator, i), 1:ndims(estimator)), 2)

        $(Symbol("rates_", f))(estimator::Estimator{<:AbstractML}) = $(Symbol("rates_", f))(estimator, 1) 

        $(Symbol("rates_", f))(estimator::Estimator{<:AbstractIndexSet}, dir::Integer) = $(Symbol("rates_", f))(estimator, (maximum(getindex.(keys(estimator), dir)) + 1) * Index(ntuple(i -> i == dir, ndims(estimator))), dir)

        function $(Symbol("rates_", f))(estimator::Estimator{<:AbstractIndexSet}, idx::Index, dir::Integer)
            m = idx[dir] - 1
            if m < 2
                return (NaN, NaN)
            else
                x = m:-1:1
                y = map(xᵢ -> $g(estimator, idx - xᵢ * Index(ntuple(i -> i == dir, ndims(estimator)))), 1:m)
                idcs = .!isnan.(y)
                θ = interp1(view(x, idcs), log2.(abs.(view(y, idcs))))
                return tuple(θ...)
            end
        end

        function $(Symbol("_regress_", g))(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index)
            p = broadcast(dir->$(Symbol("rates_", f))(estimator, index, dir), 1:ndims(estimator))
            estimates = broadcast(dir->2^(p[dir][1]+index[dir]*p[dir][2]), 1:ndims(estimator))
            estimate = mean(filter(!isnan, estimates))
            if isnan(estimate)
                p = interp($g, estimator)
                return 2^(p[1]+sum(p[2:end].*index.I))
            else
                return estimate
            end
        end
    end
end

#
# regression
#
function interp1(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    A = hcat(ones(eltype(x), length(x)), x)
    A\y
end

function interp(f::Function, estimator::Estimator)
    idx_set = filter(i -> !isempty(samples(estimator)[1][i]), CartesianIndices(size(samples(estimator)[1])))
    A = [i == 0 ? 1 : getindex(index - one(index), i) for index in idx_set, i in 0:ndims(estimator)]	
    y = map(i -> log2(f(estimator, i - one(i))), idx_set)
	try
    	return A\y
	catch e
		return fill(NaN, length(y))
	end
end

regress_mean(estimator, index) = _regress_mean(estimator, index)
regress_var(estimator, index) = _regress_var(estimator, index)
regress_cost(estimator, index) = estimator[:cost_model] isa EmptyFunction ? _regress_cost(estimator, index) : estimator[:cost_model](index)

regress_nb_of_samples(estimator::Estimator{<:SL}, index_set, ϵ::Real, θ::Real, L::Integer) = Dict(Level(0) => estimator[:nb_of_warm_up_samples])

regress_nb_of_samples(estimator::Estimator{<:AbstractIndexSet, <:QMC}, index_set, ϵ::Real, θ::Real, L::Integer) = Dict(index => estimator[:nb_of_warm_up_samples] for index in index_set)

function regress_nb_of_samples(estimator::Estimator, index_set, ϵ::Real, θ::Real, L::Integer)
    if estimator[:do_regression] && L > 2
        return _regress_nb_of_samples(estimator, index_set, ϵ, θ)
    else
        return Dict(index => estimator[:nb_of_warm_up_samples] for index in index_set)
    end
end

function _regress_nb_of_samples(estimator::Estimator{<:AbstractIndexSet, <:MC}, index_set, ϵ::Real, θ::Real)
    vars = Dict(index => regress_var(estimator, index) for index in index_set)
    costs = Dict(index => regress_cost(estimator, index) for index in index_set)
    Σ_estimate = Σ(estimator)
    for index in keys(vars)
        Σ_estimate += sqrt(vars[index] * costs[index])
    end
    Dict(index => begin
             n_opt = optimal_nb_of_samples(ϵ, θ, vars[index], costs[index], Σ_estimate)
             max(2, min(n_opt, estimator[:nb_of_warm_up_samples])) 
         end for index in index_set)
end

function compute_splitting(estimator::Estimator, ϵ::Real)
    bias_estimate = bias(estimator, max_sz(estimator))
    a = estimator[:min_splitting]
    b = estimator[:max_splitting]

    splitting = 1 - bias_estimate^2/ϵ^2
    isnan(splitting) ? a : min(b, max(a, splitting))
end

Σ(estimator::Estimator) = sum(sqrt.(map(index -> var(estimator, index) * cost(estimator, index), keys(estimator))))

optimal_nb_of_samples(estimator::Estimator, index::Index, ϵ::Real, θ::Real) = optimal_nb_of_samples(ϵ, θ, var(estimator, index), cost(estimator, index), Σ(estimator))

optimal_nb_of_samples(ϵ::Real, θ::Real, var_estimate::Real, cost_estimate::Real, Σ_estimate::Real) = ceil(Int, 1/(θ*ϵ^2) * sqrt(var_estimate/cost_estimate) * Σ_estimate)

#
# bias computation
#
boundary(estimator::Estimator, cntr::Integer) = setdiff(get_index_set(estimator, cntr), get_index_set(estimator, cntr-1))

new_index_set(estimator::Estimator, cntr::Integer) = boundary(estimator, cntr)

bias(estimator::Estimator{<:SL}) = 0.0

bias(estimator::Estimator) = bias(estimator, sz(estimator))

function bias(estimator::Estimator, sz::Integer)
    if !isempty(boundary(estimator, sz + 1) ∩ keys(estimator)) && !estimator[:robustify_bias_estimate]
        return abs(sum(broadcast(i -> mean(estimator, i), boundary(estimator, sz + 1))))
    else
        x = max(1, sz - 2):sz
        y = Float64[log2(abs(sum(broadcast(i -> mean(estimator, i), boundary(estimator, xᵢ))))) for xᵢ in x]
        p = interp1(x, y)
        return 2^(p[1]+(sz+1)*p[2])
    end
end

#
# adaptivity
#
profit(estimator::Estimator{<:AD}, index::Index) = abs(mean(estimator, index)) / sqrt(var(estimator, index) * cost(estimator, index))^estimator[:penalization]

max_level_exceeded(estimator::Estimator{<:AD}) = isempty(setdiff(get_index_set(estimator[:max_search_space], estimator[:max_index_set_param]), keys(estimator)))

function find_index_with_max_profit(estimator::Estimator{<:AD})
    indices = collect(active_set(estimator))
    profits = [profit(estimator, index) for index in indices] 
    (max_profit, idx) = findmax(profits)
    max_index = indices[idx]
    estimator[:verbose] && print_largest_profit(estimator, max_index, max_profit, indices, profits)
    return max_index
end

function new_index_set(estimator::Estimator{<:AD}, sz::Integer)
	d = ndims(estimator)
	if isempty(active_set(estimator))
		new_indices = Set{Index{d}}()
		max_index = Index(ntuple(i -> 0, d))
		add_to_active_set(estimator, max_index)
		push!(new_indices, max_index)
	else
		max_index = find_index_with_max_profit(estimator)
		add_to_old_set(estimator, max_index)
		remove_from_active_set(estimator, max_index)
		new_indices = Set{Index{d}}()
		push!(new_indices, max_index)
		for k in 1:d
			new_index = max_index + Index(ntuple(i -> i == k, d))
			if is_admissable(estimator, new_index)
				if new_index ∈ get_index_set(estimator[:max_search_space], estimator[:max_index_set_param])
					add_to_active_set(estimator, new_index)
					push!(new_indices, new_index)
				else
					warn_max_index(estimator, max_index)
					add_to_max_index_set(estimator, new_index)
				end
			end
		end
	end
	log_adaptive_index_set(estimator, max_index)
    return new_indices
end

function bias(estimator::Estimator{<:AD}, sz::Integer)
	b(indices) = abs(sum(broadcast(i -> mean(estimator, i), collect(indices))))
	if sz != max_sz(estimator)
		return b(active_set(estimator))
	else
		return min(b(active_set(estimator)), b(boundary(estimator)))
	end
end

#
# QMC related functions
#
qoi_with_max_var(estimator::Estimator{<:AbstractIndexSet, <:QMC}) = argmax(map(n_qoi -> sum(mean(var(samples_diff(estimator, n_qoi, n_shift, index)) for n_shift in 1:estimator[:nb_of_shifts](index)) for index in keys(estimator)), 1:estimator[:nb_of_qoi]))

function next_number_of_samples(estimator, index)
    if estimator[:sample_mul_factor] == 2
        Dict(index => nextpow2(nb_of_samples(estimator, index) + 1))
    elseif estimator[:sample_mul_factor] ≤ 1
        Dict(index => nb_of_samples(estimator, index) + 1)
    else
        Dict(index => ceil(Int, nb_of_samples(estimator, index) * estimator[:sample_mul_factor]))
    end
end

function find_index_with_max_var_over_cost(estimator::Estimator{<:AbstractIndexSet, <:QMC})
    indices = collect(keys(estimator))
    vars = [varest(estimator, index) / cost(estimator, index) for index in indices]
    (max_var, idx) = findmax(vars)
    indices[idx]
end

varest(estimator::Estimator{<:AbstractIndexSet, <:QMC}, index::Index) = var(mean(samples_diff(estimator, qoi_with_max_var(estimator), n_shift, index)) for n_shift in 1:estimator[:nb_of_shifts](index), corrected=true) / estimator[:nb_of_shifts](index) 

#
# Unbiased estimation
#
qoi_with_max_var(estimator::Estimator{<:U, <:MC}) = argmax(map(n_qoi -> var(accumulator(estimator, n_qoi)) / length(accumulator(estimator, n_qoi)), 1:estimator[:nb_of_qoi]))

qoi_with_max_var(estimator::Estimator{<:U, <:QMC}) = argmax(map(n_qoi -> var(mean(accumulator(estimator, n_qoi, n_shift)) for n_shift in 1:size(accumulator(estimator), 2)) / estimator[:nb_of_qoi], 1:estimator[:nb_of_qoi]))

function next_number_of_samples(estimator::Estimator{<:U})
    n_total = sum(values(nb_of_samples(estimator)))
    n0 = max(estimator[:nb_of_warm_up_samples], n_total)
    n = estimator[:sample_mul_factor] ≤ 1 ? 1 : ceil(Int, n0 * (estimator[:sample_mul_factor] - 1))
    sample(pmf(estimator), n)
end

function sample(pmf::Dict{Index{d}, Float64}, n::Integer) where d
    u = rand(n)
    x = Vector{keytype(pmf)}(undef, n)
    psum = 0.
    for (index, p) in pmf
        for i in 1:length(x)
            if psum ≤ u[i] ≤ psum + p
                x[i] = index
            end
        end
        psum = psum + p
    end
    Dict(index => sum(x .== Ref(index)) for index in keys(pmf))
end

Geometric(p, k) = (1 - p)^k*p

function update_pmf(estimator::Estimator{<:U})
    f(estimator, index) = sqrt(var(estimator, index)/cost(estimator, index))
    p = interp(f, estimator)
	if !any(isnan.(p)) && all(p[2:end] .< -0.5)
		for index in keys(estimator)
			if isnan(f(estimator, index)) || isinf(f(estimator, index))
				set_pmf_key(estimator, index, 2^(p[1]+sum(p[2:end].*index.I)))
			else
				set_pmf_key(estimator, index, f(estimator, index))
			end
		end
		normalize!(pmf(estimator))
	end
end

function normalize!(pmf::Dict{<:Index, Float64})
    tot_sum = sum(values(pmf))
    for (key, val) in pmf
        pmf[key] /= tot_sum
    end
end

Prob(estimator::Estimator{<:U}, index::Index) = sum(val for (key, val) in pmf(estimator) if key ≥ index)

Prob(estimator::Estimator{<:U}) = Dict(index => Prob(estimator, index) for index in keys(estimator))

varest(estimator::Estimator{<:U, <:MC}) = var(accumulator(estimator, qoi_with_max_var(estimator))) / length(accumulator(estimator, qoi_with_max_var(estimator)))

varest(estimator::Estimator{<:U, <:QMC}) = var(mean(accumulator(estimator, qoi_with_max_var(estimator), n_shift)) for n_shift in 1:size(accumulator(estimator), 2)) / estimator[:nb_of_qoi]

bias(estimator::Estimator{<:U}) = 0.0

converged(estimator::Estimator{<:U}, ϵ::Real) = mse(estimator) ≤ ϵ^2

mean(estimator::Estimator{<:U}) = mean(accumulator(estimator, qoi_with_max_var(estimator)))

var(estimator::Estimator{<:U}) = var(accumulator(estimator, qoi_with_max_var(estimator)))
