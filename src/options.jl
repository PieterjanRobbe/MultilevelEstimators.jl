## options.jl : stores estimator options
#
# A type that stores estimator options shared for all Estimators.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for
# Multilevel Monte Carlo Methods (c) Pieterjan Robbe, 2019

# valid options
default_options(::AbstractIndexSet, ::AbstractSampleMethod) = [:nb_of_warm_up_samples,
                                                               :nb_of_qoi,
                                                               :qoi_with_max_var,
                                                               :nb_of_tols,
                                                               :continuation_mul_factor,
                                                               :continuate,
                                                               :save,
                                                               :save_samples,
                                                               :verbose,
                                                               :folder,
                                                               :name,
                                                               :cost_model,
                                                               :nb_of_workers,
                                                               :nb_of_uncertainties,
                                                               :max_index_set_param]
get_valid_options(index_set, sample_method) = push!(default_options(index_set, sample_method), get_valid_options(index_set)..., get_valid_options(sample_method)...)

# additional valid options for IndexSet
default_options(::AbstractIndexSet) = [:min_splitting,
                                       :max_splitting,
                                       :robustify_bias_estimate,
                                       :do_mse_splitting,
                                       :do_regression]
get_valid_options(index_set::AbstractIndexSet) = default_options(index_set)
get_valid_options(::SL) = Symbol[]
get_valid_options(index_set::AD) = push!(default_options(index_set), :max_search_space, :penalization, :acceptance_rate)
get_valid_options(::U) = [:max_search_space, :sample_mul_factor]

# additional valid options for SampleMethod
get_valid_options(::AbstractSampleMethod) = Vector{Symbol}(undef, 0)
get_valid_options(::QMC) = [:nb_of_shifts,
                            :point_generator,
                            :sample_mul_factor]

# easy access to options fields
getindex(estimator::Estimator, s::Symbol) = estimator.options[s]
