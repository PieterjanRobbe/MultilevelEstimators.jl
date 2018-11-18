## init.jl : initialize lognormal diffusion problem
#
# Returns an estimator for the 2d lognormal diffusion problem.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## init_lognormal ##
function init_lognormal(index_set::AbstractIndexSet, sample_method::AbstractSampleMethod; kwargs...)

    # read optional arguments
    args = Dict{Symbol,Any}(kwargs)

    # compute Gaussian random fields
    cov_fun = CovarianceFunction(2, get_arg(args, :covariance_function))
    indices = get_index_set(index_set, get_arg(args, :max_index_set_param)) # TODO for AD?
    m0 = get_arg(args, :nb_of_coarse_dofs)
    p = get_arg(args, :minpadding)
    grf_generator = get_arg(args, :grf_generator)
    grfs = Dict(index => compute_grf(cov_fun, grf_generator, m0, index, p) for index in indices)

    # sample function
    sample_function = (index, x) -> sample_lognormal(index, x, grfs[index], get_arg(args, :damping))

    # distributions
    s = min(3600, maximum(randdim.(collect(values(grfs)))))
    distributions = [Normal() for i in 1:s]

    # estimator
    filter_keys!(args)
    Estimator(index_set, sample_method, sample_function, distributions; args...)

end

## compute_grf ##
grid_size(index::Index) = 2 .^index
grid_size(level::Level) = (2^level[1], 2^level[1])

function compute_grf(cov_fun::CovarianceFunction, grf_generator::GaussianRandomFieldGenerator, m0, index::Index, p)
    n = m0 .* grid_size(index)
    pts = broadcast(i -> range(0, stop=1, length=i+1), n)
    compute_grf(cov_fun, grf_generator, pts, p)
end

compute_grf(cov_fun, grf_generator::GaussianRandomFieldGenerator, pts, p) = GaussianRandomField(cov_fun, grf_generator, pts...)

compute_grf(cov_fun, grf_generator::CirculantEmbedding, pts, p) = GaussianRandomField(cov_fun, grf_generator, pts..., minpadding=p)

## get_arg ##
macro get_arg(key_name, default_value)
    return eval(:(get_arg(args::Dict{Symbol, Any}, ::Val{$key_name}) = haskey(args, $key_name) ? args[$key_name] : $default_value))
end

get_arg(args::Dict{Symbol,Any}, arg::Symbol) = get_arg(args, Val(arg))

get_arg(args::Dict{Symbol,Any}, arg::Val{T}) where T = throw(ArgumentError(string("in init_lognormal, invalid key ", T, " found")))

@get_arg :nb_of_coarse_dofs 4

@get_arg :covariance_function Matern(get_arg(args, :length_scale), get_arg(args, :smoothness))

@get_arg :length_scale 0.1

@get_arg :smoothness 1

@get_arg :max_index_set_param 6

@get_arg :grf_generator CirculantEmbedding()

@get_arg :minpadding 0

@get_arg :damping 0.8

# make sure all keys in args are valid keys for Estimator
filter_keys!(args::Dict{Symbol, Any}) = isempty(args) || delete!.(Ref(args), [:nb_of_coarse_dofs, :covariance_function, :length_scale, :smoothness, :max_index_set_param, :grf_generator, :minpadding])
