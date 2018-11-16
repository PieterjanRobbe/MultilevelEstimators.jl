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

	# compute Gaussian random fields
	cov_fun = CovarianceFunction(2, get_arg(args, :covariance_function))

end

## get_arg ##
function make_get_arg(key_name, default_value)
	@show key_name
	#@show eval(:($key_name))
	key = :length_scale
	@show key
	key_name = Expr(:quote, key_name)
	@show key_name
	#sym2 = eval(sym)
	#@show sym2
	#@show eval(:(:($$(key))))
	eval(
		 quote
			 #get_arg(args::Dict{Symbol, Any}, ::Val{$(key_name)}) = haskey(args, $(key_name)) ? args[$(key_name)] : $(default_value)
			 get_arg(args::Dict{Symbol, Any}, ::Val{$(key_name)}) = typeof($key_name)
		 end)
end


get_arg(args::Dict{Symbol,Any}, arg::Symbol) = get_arg(args, Val(arg))

get_arg(args::Dict{Symbol,Any}, arg::Val{T}) where T = throw(ArgumentError(string("in init_lognormal, invalid key ", T, " found")))

get_arg(args::Dict{Symbol, Any}, ::Val{:nb_of_coarse_dof}) = haskey(args, :nb_of_coarse_dof) ? args[:nb_of_coarse_dof] : 2

get_arg(args::Dict{Symbol, Any}, ::Val{:covariance_function}) = haskey(args, :covariance_function) ? args[:covariance_function] : Matern(get_arg(args, :length_scale), get_arg(args, :smoothness))


make_get_arg(:length_scale, 1)

#get_arg(args::Dict{Symbol, Any}, ::Val{:length_scale}) = haskey(args, :length_scale) ? args[:length_scale] : 1

get_arg(args::Dict{Symbol, Any}, ::Val{:smoothness}) = haskey(args, :smoothness) ? args[:smoothness] : 1
