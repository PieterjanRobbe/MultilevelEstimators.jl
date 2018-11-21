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
    for key in valid_keys(index_set, sample_method)
        parse!(index_set, sample_method, settings, key)
    end

    # create options and internals
    options = EstimatorOptions(settings)
    internals = EstimatorInternals(index_set, sample_method, settings)

    Estimator{typeof(index_set), typeof(sample_method), typeof(distributions), typeof(options), typeof(internals)}(index_set, sample_function, distributions, options, internals)
end

Estimator(index_set::AbstractIndexSet, sample_method::AbstractSampleMethod, sample_function::Function, distribution::AbstractDistribution; kwargs...) = Estimator(index_set, sample_method, sample_function, [distribution]; kwargs...)

## check valid keys ##
function check_valid_keys(settings::Dict, index_set::AbstractIndexSet, sample_method::AbstractSampleMethod)
    for key in keys(settings)
        key ∉ valid_keys(index_set, sample_method) && throw(ArgumentError(string("in Estimator, invalid key ", key, " found")))
    end
end

## valid keys ##
valid_keys(index_set::AbstractIndexSet, sample_method::AbstractSampleMethod) = vcat(valid_keys(), valid_keys(index_set), valid_keys(sample_method))

valid_keys() = collect(fieldnames(EstimatorOptions))

valid_keys(::AbstractSampleMethod) = Symbol[]
valid_keys(::QMC) = [:nb_of_shifts, :point_generator, :sample_mul_factor]
valid_keys(::MC) = [:do_regression]

valid_keys(::AbstractIndexSet) = Symbol[]
valid_keys(::AbstractAD) = [:max_search_space]
valid_keys(::MG) = [:sample_mul_factor]
valid_keys(::MG{<:AD}) = [:sample_mul_factor, :max_search_space]

## print methods ##
show(io::IO, estimator::Estimator{I,S}) where {I<:AbstractIndexSet, S<:AbstractSampleMethod} = print(io, string("Estimator{", estimator.index_set, ", ", S, "}"))

## getters and setters ##
push!(estimator::Estimator, index::Index) = push!(estimator.internals.current_index_set, index)

nb_of_warm_up_samples(estimator::Estimator) = estimator.options.nb_of_warm_up_samples

contains_samples_at_index(estimator::Estimator, index::Index) = isassigned(estimator.internals.samples_diff, 1) && haskey(estimator.internals.samples_diff[1], index)

total_work(estimator::Estimator) = estimator.internals.total_work
total_work(estimator::Estimator, index::Index) = get(total_work(estimator), index, nothing)

nb_of_samples(estimator::Estimator) = estimator.internals.nb_of_samples
nb_of_samples(estimator::Estimator, index::Index) = get(nb_of_samples(estimator), index, 0)

cost_model(estimator::Estimator) = estimator.internals.cost_model
cost_model(estimator::Estimator, index::Index) = cost_model(estimator)(index)

nb_of_workers(estimator::Estimator, index::Index) = estimator.options.nb_of_workers(index)

nb_of_shifts(estimator::Estimator{<:AbstractIndexSet, <:QMC}, index::Index) = estimator.internals.nb_of_shifts(index)

nb_of_qoi(estimator::Estimator) = estimator.options.nb_of_qoi

name(estimator) = estimator.options.name # TODO automate?

shortname(estimator) = first(split(estimator.options.name, "."))

function add_index(estimator::Estimator, index::Index)
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

distributions(estimator::Estimator) = estimator.distributions

stochastic_dim(estimator::Estimator) = length(distributions(estimator::Estimator))

nb_of_tols(estimator::Estimator) = estimator.options.nb_of_tols

continuation_mul_factor(estimator::Estimator) = estimator.options.continuation_mul_factor

continuate(estimator::Estimator) = estimator.options.continuate

verbose(estimator::Estimator) = estimator.options.verbose

get_tols(estimator::Estimator, tol::T) where T<:Real = continuate(estimator) ? continuation_mul_factor(estimator).^(nb_of_tols(estimator)-1:-1:0)*tol : T[tol] 

samples(estimator::Estimator) = estimator.internals.samples
samples(estimator::Estimator{<:AbstractIndexSet, <:MC}, n_qoi::Integer) = getindex(samples(estimator), n_qoi)
samples(estimator::Estimator{<:AbstractIndexSet, <:MC}, n_qoi::Integer, index::Index) = get(samples(estimator, n_qoi), index, nothing)

samples_diff(estimator::Estimator) = estimator.internals.samples_diff
samples_diff(estimator::Estimator{<:AbstractIndexSet, <:MC}, n_qoi::Integer) = getindex(samples_diff(estimator), n_qoi)
samples_diff(estimator::Estimator{<:AbstractIndexSet, <:MC}, n_qoi::Integer, index::Index) = get(samples_diff(estimator, n_qoi), index, nothing)

append_samples!(estimator::Estimator{<:AbstractIndexSet, <:MC}, n_qoi::Integer, index::Index, samples_to_append) = append!(estimator.internals.samples[n_qoi][index], samples_to_append)
append_samples_diff!(estimator::Estimator{<:AbstractIndexSet, <:MC}, n_qoi::Integer, index::Index, samples_to_append) = append!(estimator.internals.samples_diff[n_qoi][index], samples_to_append)

update_nb_of_samples!(estimator::Estimator, index::Index, nb_of_samples::Integer) = estimator.internals.nb_of_samples[index] += nb_of_samples

update_total_work!(estimator::Estimator, index::Index, time::Real) = total_work(estimator)[index] += cost_model(estimator) isa EmptyFunction ? time : estimator.internals.cost_model(index)

keys(estimator::Estimator) = sort(collect(estimator.internals.current_index_set))
