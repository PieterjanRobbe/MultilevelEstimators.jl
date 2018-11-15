## estimators.jl : main Estimator type
#
# Representation of the main Estimator type.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## Estimator ##
struct Estimator{I<:AbstractIndexSet, S<:AbstractSampleMethod,D}
    index_set::I
    sample_function::Function
    distributions::D
end

function Estimator(index_set::AbstractIndexSet, sample_method::AbstractSampleMethod, sample_function::Function, distributions::Vector{<:AbstractDistribution}; kwargs...)

    # read optional arguments
    settings = Dict{Symbol,Any}(kwargs)
    check_valid_keys(settings, index_set, sample_method)
    
    # default settings
    for key in valid_keys(index_set, sample_method)
        parse!(index_set, sample_method, settings, key)
    end

    # create options, base_estimator...
    #display(settings)
    
    Estimator{typeof(index_set), typeof(sample_method), typeof(distributions)}(index_set, sample_function, distributions)
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
valid_keys(::AbstractIndexSet) = Symbol[]
valid_keys(::AbstractSampleMethod) = Symbol[]

# TODO replace this by fieldnames(EstimatorOptions) ?
valid_keys() = [:nb_of_warm_up_samples, :nb_of_qoi, :continuate, :nb_of_tols, :continuation_mul_factor, :folder, :name, :save_samples, :user_data, :verbose, :cost_model, :robustify_bias_estimate, :do_mse_splitting, :sample_mul_factor, :max_index_set_param, :nb_of_workers]
valid_keys(::QMC) = [:nb_of_shifts, :point_generator, :sample_mul_factor]
valid_keys(::MC) = [:do_regression]
valid_keys(::AbstractAD) = [:max_search_space]
valid_keys(::MG) = [:sample_mul_factor]
valid_keys(::MG{<:AD}) = [:sample_mul_factor, :max_search_space]

## print methods ##
show(io::IO, estimator::Estimator{I,S}) where {I<:AbstractIndexSet, S<:AbstractSampleMethod} = print(io, string("Estimator{", estimator.index_set, ", ", S, "}"))
