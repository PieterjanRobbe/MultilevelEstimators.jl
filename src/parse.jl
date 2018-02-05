## parse.jl : parse all inputs and set the corresponding sampler options

"""
    Sampler(kwargs...)

Create a multilevel `Sampler` object using the options specified in `kwargs`.

# Keywords
method
number_generator
sample_function

user_data
cost_model
use_sample_run_time
nb_of_qoi
nb_of_warm_up_samples
nb_of_tol
maximum_level
maximum_index
valid_indices
low_mem
"""
function Sampler(kwargs...)
    settings = Dict(kwargs)

	# TODO : either provide sample_function or ml_sample_function!!!
	# provide default ml_sample_function for MLMC, MIMC that returns an approx on coarse and fine scale
	# and uses sample_fun

    # required
    for key in [:method, :number_generator, :sample_function]
        haskey(settings,key) || throw(ArgumentError("required key $(key) not provided!"))
    end
    
    # create basic sampler
    sampler = Sampler(settings[:method],settings[:number_genertor],settings[:sample_function],settings[:user_data]) 

    # optional
    for (key, value) in settings
        set_key!(sampler, key, value)
    end
end

set_key!(::Sampler, ::Type{Val{T}}, ::Any) where T = throw(ArgumentError("could not parse key $(key)!"))

## REQUIRED ##

set_key!(sampler,::Type{Val{:method}},method::IndexSet) = nothing
set_key!(sampler,::Type{Val{:number_generator}},number_generator::NumberGenerator) = nothing
set_key!(sampler,::Type{Val{:sample_function}},sample_function::Function) = nothing

set_key!(sampler,::Type{Val{:user_data}},user_data) = nothing

## OPTIONAL ##

# cost_model
function set_key!(sampler,::Type{Val{:cost_model}},cost_model::Function)
    sampler.use_sample_run_time == true && warn("conflicting options use_sample_run_time and cost_model detected, using cost_model")
    sampler.use_sample_run_time = false
    sampler.cost_model = cost_model
end

# use_sample_run_time
function set_key!(sampler,::Type{Val{:use_sample_run_time}},use_sample_run_time::Bool)
    ( use_sample_run_time == && typeof(sampler.cost_model) <: EmptyCostModel ) || warn("conflicting options use_sample_run_time and cost_model detected, using use_sample_run_time")
    sampler.use_sample_run_time = use_sample_run_time
    use_sample_run_time && sampler.cost_model = EmptyCostModel()
end

# nb_of_qoi
function set_key!(sampler,::Type{Val{:nb_of_qoi}},nb_of_qoi::N where N<:Integer)
    nb_of_qoi > 0 || throw(ArgumentError("nb_of_qoi must be positive, got $(nb_of_qoi)"))
    sampler.nb_of_qoi = nb_of_qoi
end

# nb_of_warm_up_samples
function set_key!(sampler,::Type{Val{:nb_of_warm_up_samples}},nb_of_warm_up_samples::N where N<:Integer)
    nb_of_warm_up_samples > 0 || throw(ArgumentError("nb_of_warm_up_samples must be positive, got $(nb_of_warm_up_samples)!"))
    sampler.nb_of_warm_up_samples = nb_of_warm_up_samples
end

# nb_of_tol
function set_key!(sampler,::Type{Val{:nb_of_tol}},nb_of_tol::N where N<:Integer)
    nb_of_tol > 0 || throw(ArgumentError("nb_of_tol must be positive, got $(nb_of_tol)!"))
    sampler.nb_of_tol = nb_of_tol
end

# maximum_level
function set_key!(sampler,::Type{Val{:maximum_level}},maximum_level::Level)
    typeof(sampler.method) <: ML || throw(ArgumentError("cannot set maximum level if sampler uses a $(sampler.method)"))
    isempty(sampler.valid_indices) || warn("conflicting options maximum_level and valid_indices detected, using maximum_level")
    sampler.valid_indices = get_indexset(ML(),maximum_level[1])
end

# maximum_index
function set_key!(sampler,::Type{Val{:maximum_level}},maximum_index::Index{d})
    typeof(sampler.method) <: Union{SL,ML} && throw(ArgumentError("cannot set maximum index if sampler uses a $(sampler.method)"))
    ndims(sampler.method) == d || throw(DimensionMismatch("number of dimensions in maximum_index does not agree with number of dimensions of the sampler (expected $(ndims(sampler.method)), got $(d))"))
    isempty(sampler.valid_indices) || warn("conflicting options maximum_index and valid_indices detected, using maximum_index")
    N = maximum(maximum_index)
    idxset = FT(ndims(sampler.method),δ=1./maximum_index*N)
    sampler.valid_indices = get_indexset(idxset,N)
end

# valid_indices
function set_key!(sampler,::Type{Val{:valid_indices}},valid_indices::Set{Index{d}}) where {d}
    ndims(sampler.method) == d || throw(DimensionMismatch("number of dimensions in valid_indices does not agree with number of dimensions of the sampler (expected $(ndims(sampler.method)), got $(d))"))
    isempty(sampler.valid_indices) || warn("setting valid_indices, but already specified")
    is_valid_index_set(valid_indices) || throw(ArgumentError("valid_indices is not a valid index set"))
    sampler.valid_indices = valid_indices
end

# low_mem
function set_key!(sampler,::Type{Val{:low_mem}},low_mem::Bool)
    sampler.low_mem = low_mem
end

## nullables
type EmptyCostModel <: Function end








#####################################
#####################################
#####################################
#####################################
# TODO define some standard cost models: \gamma and d, and provide multi_index_cost()
# TODO define use_sample_run_time = true (default), unless a cost_model is provided
# TODO default value is Nstar = convert(Int64,ceil(32/nshifts(numberGenerators[Index(zeros(Int64,d))])))
# TODO avoid the "safety" keyword, only use doubling algorithm in QMC setting
# TODO avoid the continuate keyword, see algorithm implementation AND use the nb_of_tol keyword
# TODO how about "real" continuation ?? specify k0 and k1 (trust parameters)
# TODO what about the old store_samples_0 key

#=
DEFAULTS
    generatorStates = Dict{Index{d,Vector{Int64}},Int64}()

    samples = Dict{Index{d,Vector{Int64}},Array{Float64,3}}()
    samples0 = Dict{Index{d,Vector{Int64}},Array{Float64,3}}()

    T = Dict{Index{d,Vector{Int64}},Float64}()
    E = Dict{Index{d,Vector{Int64}},Vector{Float64}}()
    V = Dict{Index{d,Vector{Int64}},Vector{Float64}}()
    Vf = Dict{Index{d,Vector{Int64}},Vector{Float64}}()
    Wst = Dict{Index{d,Vector{Int64}},Float64}()
    W = Dict{Index{d,Vector{Int64}},Float64}()
    Vest = Dict{Index{d,Vector{Int64}},Vector{Float64}}()
    P = Dict{Index{d,Vector{Int64}},Float64}()

    return Sampler{d,typeof(indexSet),typeof(numberGenerators),typeof(gaussianFieldSampler),typeof(generatorStates),typeof(samples),typeof(T),typeof(E),typeof(userType),typeof(max_indexset)}(
                                                                                                                                                                                               indexSet, numberGenerators, sampleFunction, mmaxL, costModel, Z, Nstar, gaussianFieldSampler, useTime,
                                                                                                                                                                                               safety, continuate, nTOL, k, showInfo, ioStream, storeSamples0, procMap, userType, max_indexset, 
                                                                                                                                                                                               generatorStates, samples, samples0, T,E,V,Vf,Wst,W,Vest,P, ml_sample_fun,reuseSamples, Dict{Index{d,Vector{Int64}},Int64}(),is_cauchy_schwarz, use_batches, giles_multigrid)
end

=# 
