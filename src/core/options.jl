## options.jl : stores estimator options
#
# A type that stores estimator options shared for all Estimators.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## EstimatorOptions ##
struct EstimatorOptions{N<:Integer, T<:Real}
    nb_of_warm_up_samples::N
    nb_of_qoi::N
    max_index_set_param::N
    nb_of_tols::N

    continuation_mul_factor::T
    splitting::T

    continuate::Bool
    save_samples::Bool
    robustify_bias_estimate::Bool
    do_mse_splitting::Bool
    verbose::Bool

    folder::String
    name::String

    cost_model::Function
    nb_of_workers::Function
end

EstimatorOptions(settings::Dict{Symbol,Any}) = 
EstimatorOptions(
                 promote([settings[name] for name in fieldnames(EstimatorOptions)[1:4]]...)..., # promotion of N
                 promote([settings[name] for name in fieldnames(EstimatorOptions)[5:5]]...)..., # promotion of T
                 [settings[name] for name in fieldnames(EstimatorOptions)[6:end]]... # other fields
                )
