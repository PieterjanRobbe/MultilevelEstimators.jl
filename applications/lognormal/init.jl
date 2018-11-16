## init.jl : initialize lognormal diffusion problem
#
# Returns an estimator for the 2d lognormal diffusion problem.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## init_lognormal ##
function init_lognormal(index_set::MultilevelEstimators.AbstractIndexSet, sample_method::MultilevelEstimators.AbstractSampleMethod; kwargs...)

    # read optional arguments
    args = Dict{Symbol,Any}(kwargs)

    # level/index hierarchy
    @show get_arg(args, :length_scale)
    #@show get_arg(args, :covariance_function)

    # compute Gaussian random fields
    #cov_fun = CovarianceFunction(2, get_arg(args, :covariance_function))

end

## get_arg ##
macro get_arg(key_name, default_value)
    return eval(:(get_arg(args::Dict{Symbol, Any}, ::Val{$key_name}) = haskey(args, $key_name) ? args[$key_name] : $default_value))
end


get_arg(args::Dict{Symbol,Any}, arg::Symbol) = get_arg(args, Val(arg))

get_arg(args::Dict{Symbol,Any}, arg::Val{T}) where T = throw(ArgumentError(string("in init_lognormal, invalid key ", T, " found")))

@get_arg :nb_of_coarse_dofs 2

@get_arg :covariance_function Matern(get_arg(args, :length_scale), get_arg(args, :smoothness))

@get_arg :length_scale 1

@get_arg :smoothness 1
