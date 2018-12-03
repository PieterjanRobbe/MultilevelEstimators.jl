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
    args[:index_set] = index_set

    # compute Gaussian random fields
    cov_fun = CovarianceFunction(2, get_arg(args, :covariance_function))
    indices = get_max_index_set(index_set, args)
    m0 = get_arg(args, :nb_of_coarse_dofs)
    p = get_arg(args, :minpadding)
    grf_generator = get_arg(args, :grf_generator)
    grfs = Dict(index => compute_grf(cov_fun, grf_generator, m0, index, p) for index in indices)

    # sample function
    damping = get_arg(args, :damping)
    qoi = get_arg(args, :qoi)
    solver = get_arg(args, :solver)
    reuse = get_arg(args, :reuse) ? Reuse() : NoReuse()
    sample_function = (index, x) -> sample_lognormal(index, x, grfs[index], damping, qoi, solver, reuse)

    # distributions
    s = maximum(randdim.(collect(values(grfs))))
    distributions = [Normal() for i in 1:s]

    # set nb_of_uncertainties for more efficient sampling
    args[:nb_of_uncertainties] = index -> randdim(grfs[index])

    # set max_index_set_param
    args[:max_index_set_param] = get_arg(args, :max_index_set_param)

    # estimator
    filter_keys!(args)
    Estimator(index_set, sample_method, sample_function, distributions; args...)

end

## get_max_index_set ##
get_max_index_set(index_set::AbstractIndexSet, args::Dict{Symbol, Any}) = get_index_set(index_set, get_arg(args, :max_index_set_param))
get_max_index_set(::SL, args::Dict{Symbol, Any}) = [Level(0)]
# TODO for AD ?

## compute_grf ##
grid_size(index::Index) = 2 .^index
grid_size(level::Level) = (2^level[1], 2^level[1])

function compute_grf(cov_fun::CovarianceFunction, grf_generator::GaussianRandomFieldGenerator, m0, index::Index, p)
    n = m0 .* grid_size(index)
    pts = broadcast(i -> range(0, stop=1, length=i+1), n)
    compute_grf(cov_fun, grf_generator, pts, p)
end

compute_grf(cov_fun, grf_generator::GaussianRandomFieldGenerator, pts, p) = GaussianRandomField(cov_fun, grf_generator, pts...)

compute_grf(cov_fun, grf_generator::CirculantEmbedding, pts, p) = GaussianRandomField(cov_fun, grf_generator, pts..., minpadding=p, measure=false)

## get_arg ##
macro get_arg(key_name, default_value)
    return eval(:(get_arg(args::Dict{Symbol, Any}, ::Val{$key_name}) = haskey(args, $key_name) ? args[$key_name] : $default_value))
end

get_arg(args::Dict{Symbol,Any}, arg::Symbol) = get_arg(args, Val(arg))

get_arg(args::Dict{Symbol,Any}, arg::Val{T}) where T = throw(ArgumentError(string("in init_lognormal, invalid key ", T, " found")))

@get_arg :nb_of_coarse_dofs args[:index_set] isa SL ? 4*2^get_arg(args, :max_index_set_param) : 4

@get_arg :covariance_function Matern(get_arg(args, :length_scale), get_arg(args, :smoothness))

@get_arg :length_scale 0.1

@get_arg :smoothness 1

@get_arg :max_index_set_param 6

@get_arg :grf_generator CirculantEmbedding()

@get_arg :minpadding 0

@get_arg :damping 0.8

abstract type AbstractQoi end

struct Qoi1 <:AbstractQoi end

struct Qoi2 <:AbstractQoi end

struct Qoi3 <:AbstractQoi end

@get_arg :qoi Qoi1()

abstract type AbstractSolver end

struct MGSolver <: AbstractSolver end

struct MSGSolver <: AbstractSolver end

@get_arg :solver args[:index_set] isa MG{<:MultilevelEstimators.MI} ? MSGSolver() : MGSolver()

abstract type AbstractReuse end

struct Reuse <: AbstractReuse end

struct NoReuse <: AbstractReuse end

@get_arg :reuse args[:index_set] isa MG ? true : false

# make sure all keys in args are valid keys for Estimator
filter_keys!(args::Dict{Symbol, Any}) = isempty(args) || delete!.(Ref(args), [:nb_of_coarse_dofs, :covariance_function, :length_scale, :smoothness, :grf_generator, :minpadding, :index_set, :qoi, :damping, :solver, :reuse])
