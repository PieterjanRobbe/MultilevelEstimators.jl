## estimators.jl : main Estimator type
#
# Representation of the main Estimator type.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## Estimator ##
"""
```julia
Estimator(index_set, sample_method, sample_function, distributions; kwargs...)
```
Create an `Estimator` with index set `index_set` and sample method `sample_method` for the expected value of the quantity of interest returned by `sample_function`, and where `distributions` is the uncertainty on the input parameters.

# Examples
```jldoctest
julia>

```
"""
struct Estimator{I<:AbstractIndexSet, S<:AbstractSampleMethod, D, O, N}
    index_set::I
    sample_function::Function
    distributions::D
    options::O
    internals::N
end

function Estimator(index_set::AbstractIndexSet, sample_method::AbstractSampleMethod, sample_function::Function, distributions::AbstractVector{<:AbstractDistribution}; kwargs...)

    # read optional arguments
    settings = Dict{Symbol,Any}(kwargs)
    check_valid_keys(settings, index_set, sample_method)

    # default settings
	settings[:distributions] = distributions
	for key in fieldnames(EstimatorOptions)
        parse!(index_set, sample_method, settings, key)
    end
	parse_additional_options!(index_set, sample_method, settings)

    # create options and internals
    options = EstimatorOptions(settings)
    internals = EstimatorInternals(index_set, sample_method, settings)

    Estimator{typeof(index_set), typeof(sample_method), typeof(distributions), typeof(options), typeof(internals)}(index_set, sample_function, distributions, options, internals)
end

# for a single uncertain parameter
Estimator(index_set::AbstractIndexSet, sample_method::AbstractSampleMethod, sample_function::Function, distribution::AbstractDistribution; kwargs...) = Estimator(index_set, sample_method, sample_function, [distribution]; kwargs...)

## check valid keys ##
function check_valid_keys(settings::Dict, index_set::AbstractIndexSet, sample_method::AbstractSampleMethod)
    for key in keys(settings)
        key ∉ valid_keys(index_set, sample_method) && throw(ArgumentError(string("in Estimator, invalid key ", key, " found")))
    end
end

## parse additional options ##
function parse_additional_options!(index_set::T1, sample_method::T2, settings::Dict{Symbol, Any}) where {T1<:AbstractIndexSet, T2<:AbstractSampleMethod}
	if T2 <: QMC 
		parse!(index_set, sample_method, settings, :point_generator)
	end
	if T1 <: AbstractAD
		parse!(index_set, sample_method, settings, :max_search_space)
	end
end

## valid keys ##
function valid_keys(::I, ::S) where {I<:AbstractIndexSet, S<:AbstractSampleMethod}
	v = [:nb_of_warm_up_samples,
		 :nb_of_qoi,
		 :nb_of_tols,
		 :continuation_mul_factor,
		 :continuate,
		 :save_samples,
		 :verbose,
		 :folder,
		 :name,
		 :cost_model,
		 :nb_of_workers]
	I <: Union{ML, MI} && push!(v,
								:max_index_set_param,
								:min_splitting,
								:max_splitting,
								:robustify_bias_estimate,
								:do_mse_splitting,
								:do_regression)
	I <:AbstractAD && push!(v, :max_search_space)
	I <:MG && push!(v, :sample_mul_factor)
	S <: QMC && push!(v, :nb_of_shifts,
					  :point_generator,
					  :sample_mul_factor)
	return v
end

## print methods ##
show(io::IO, estimator::Estimator{I,S}) where {I<:AbstractIndexSet, S<:AbstractSampleMethod} = print(io, string("Estimator{", estimator.index_set, ", ", S, "}"))

## getters for options ##
for f in fieldnames(EstimatorOptions)
    eval(
         quote
             $f(estimator::Estimator) = estimator.options.$f
         end)
end

## getters and setters for internals ##
contains_samples_at_index(estimator::Estimator, index::Index) = isassigned(estimator.internals.samples_diff, 1) && haskey(estimator.internals.samples_diff[1], index)

total_work(estimator::Estimator) = estimator.internals.total_work
total_work(estimator::Estimator, index::Index) = get(total_work(estimator), index, NaN)
update_total_work!(estimator::Estimator, index::Index, time::Real, n::Integer) = total_work(estimator)[index] += cost_model(estimator) isa EmptyFunction ? time : n * cost_model(estimator)(index)

nb_of_samples(estimator::Estimator) = estimator.internals.nb_of_samples
nb_of_samples(estimator::Estimator, index::Index) = get(nb_of_samples(estimator), index, 0)
update_nb_of_samples!(estimator::Estimator, index::Index, nb_of_samples::Integer) = estimator.internals.nb_of_samples[index] += nb_of_samples

cost_model(estimator::Estimator, index::Index) = cost_model(estimator)(index)

nb_of_workers(estimator::Estimator, index::Index) = nb_of_workers(estimator)(index)

nb_of_shifts(estimator::Estimator{<:AbstractIndexSet, <:QMC}, index::Index) = estimator.internals.nb_of_shifts(index)

shortname(estimator) = first(split(estimator.options.name, "."))

distributions(estimator::Estimator) = estimator.distributions

stochastic_dim(estimator::Estimator) = length(distributions(estimator::Estimator))

get_tols(estimator::Estimator, tol::T) where T<:Real = continuate(estimator) ? continuation_mul_factor(estimator).^(nb_of_tols(estimator)-1:-1:0)*tol : T[tol] 

samples(estimator::Estimator) = estimator.internals.samples
samples(estimator::Estimator{<:AbstractIndexSet, <:MC}, n_qoi::Integer) = getindex(samples(estimator), n_qoi)
samples(estimator::Estimator{<:AbstractIndexSet, <:MC}, n_qoi::Integer, index::Index) = get(samples(estimator, n_qoi), index, NaN)

samples_diff(estimator::Estimator) = estimator.internals.samples_diff
samples_diff(estimator::Estimator{<:AbstractIndexSet, <:MC}, n_qoi::Integer) = getindex(samples_diff(estimator), n_qoi)
samples_diff(estimator::Estimator{<:AbstractIndexSet, <:MC}, n_qoi::Integer, index::Index) = get(samples_diff(estimator, n_qoi), index, NaN)

append_samples!(estimator::Estimator{<:AbstractIndexSet, <:MC}, n_qoi::Integer, index::Index, samples_to_append) = append!(estimator.internals.samples[n_qoi][index], samples_to_append)
append_samples_diff!(estimator::Estimator{<:AbstractIndexSet, <:MC}, n_qoi::Integer, index::Index, samples_to_append) = append!(estimator.internals.samples_diff[n_qoi][index], samples_to_append)

keys(estimator::Estimator) = sort(collect(current_index_set(estimator)))
all_keys(estimator::Estimator) = sort(collect(keys(first(samples(estimator)))))

current_index_set(estimator::Estimator) = estimator.internals.current_index_set

ndims(::Estimator{<:AbstractIndexSet{d}}) where d = d 

get_index_set(estimator::Estimator, level) = get_index_set(estimator.index_set, level)

boundary(estimator::Estimator, cntr) = cntr == 0 ? get_index_set(estimator, cntr) : setdiff(get_index_set(estimator, cntr), get_index_set(estimator, cntr-1))

sz(estimator::Estimator) = estimator.internals.index_set_size.sz
function set_sz(estimator::Estimator, n)
    estimator.internals.index_set_size.sz = n
    estimator.internals.index_set_size.max_sz = max(sz(estimator), max_sz(estimator))
end
max_sz(estimator::Estimator) = estimator.internals.index_set_size.max_sz

## clear ##
function clear(estimator::Estimator)
    for index in keys(estimator)
        delete!(current_index_set(estimator), index)
    end
end

## push index to estimator ##
push!(estimator::Estimator, index::Index) = push!(estimator.internals.current_index_set, index)

function add_index(estimator::Estimator, index::Index) # this is called in "sample"
    estimator.internals.nb_of_samples[index] = zero(valtype(nb_of_samples(estimator)))
    estimator.internals.total_work[index] = zero(valtype(total_work(estimator)))
    for S in [samples_diff(estimator), samples(estimator)]
        for q in 1:nb_of_qoi(estimator)
            if !isassigned(S, q)
                S[q] = eltype(S)()
            end
            S[q][index] = valtype(S[q])(undef, 0)
        end
    end
end
