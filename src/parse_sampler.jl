## parse_sampler.jl : parse a sampler from an input dict with simulation settings

#TODO move this to Sampler, make empty sampler type
# parse all inputs from the input dict and set the appropriate sampler settings
function Sampler(settings::Dict{Symbol,Any})
    sampler = Sampler()

    ## required ##
    for key in [:index_set, :number_generator, :sample_function]
        haskey(settings,key) || throw(ArgumentError("required key $(key) not provided!"))
        set_key(samper, key, settings[key])
    end

    ## optional ##
    for (key, value) in settings
        set_key(sampler, key, value)
    end
end

set_key(::Sampler, ::Type{Val{T}}, ::Any) where T = throw(ArgumentError("could not parse key $(key)!"))

# index_set
function set_key(sampler,::Type{Val{:index_set}},index_set::IndexSet)
    sampler.index_set = index_set
end

# number_generator
function set_key(sampler,::Type{Val{:number_generator}},number_generator::NumberGenerator)
    # TODO generate new instance of number_generator for each new index added 
    sampler.number_generator = number_generator
end

# sample_function
function set_key(sampler,::Type{Val{:sample_function}},sample_function::Function)
    sampler.sample_function = sample_function
end

# ml_sample_fun
function set_key(sampler,::Type{Val{:ml_sample_fun}},ml_sample_fun::Function)
    sampler.ml_sample_fun = ml_sample_fun
end

# maximum_level
function set_key(sampler,::Type{Val{:maximum_level}},maximum_level::N where N<:Integer)
    maximum_level ⩾ 0 || throw(ArgumentError("maximum_level must be positive, got $(maximum_level)!"))
    sampler.maximum_level = typeof(sampler.index_set) == SL ? 0 : maximum_level = 0
end

# cost_model
function set_key(sampler,::Type{Val{:cost_model}},cost_model::Function)
    # TODO define some standard cost models: \gamma and d, and provide multi_index_cost()
    # TODO define use_sample_run_time = true (default), unless a cost_model is provided
    sampler.use_sample_run_time = false
    sampler.cost_model = cost_model
end

# nb_of_qoi
function set_key(sampler,::Type{Val{:nb_of_qoi}},nb_of_qoi::N where N<:Integer)
    nb_of_qoi > 0 || throw(ArgumentError("nb_of_qoi must be positive, got $(nb_of_qoi)!"))
    sampler.nb_of_qoi = nb_of_qoi
end

# nb_of_warm_up_samples
function set_key(sampler,::Type{Val{:nb_of_warm_up_samples}},nb_of_warm_up_samples::N where N<:Integer)
    # TODO default value is Nstar = convert(Int64,ceil(32/nshifts(numberGenerators[Index(zeros(Int64,d))])))
    nb_of_warm_up_samples > 0 || throw(ArgumentError("nb_of_warm_up_samples must be positive, got $(nb_of_warm_up_samples)!"))
end

# gaussian_random_field_sampler
function set_key(sampler,::Type{Val{:gaussian_random_field_sampler}},gaussian_random_field_sampler::GaussianRandomFieldSampler)
    # TODO should somehow check that maximum_level agrees with number eigenfunctions
    # TODO and how do we handle the HC case? the AD case?
    # TODO let the GRF type hold multiple KL expansions to cope with heat exchanger example
    # TODO EmptySampler is default, rename into EmptyGaussianRoandomFieldSampler
    sampler.gaussian_random_field_sampler = gaussian_random_field_sampler
end

# TODO avoid the "safety" keyword, only use doubling algorithm in QMC setting
# TODO avoid the continuate keyword, see algorithm implementation AND use the nb_of_tol keyword

# nb_of_tol
function set_key(sampler,::Type{Val{:nb_of_tol}},nb_of_tol::N where N<:Integer)
    nb_of_tol > 0 || throw(ArgumentError("nb_of_tol must be positive, get $(nb_of_tol)!"))
    sampler.nb_of_tol = nb_of_tol
end
 
# TODO how about "real" continuation ?? specify k0 and k1 (trust parameters)

# store_samples
function set_key(sampler,::Type{Val{:store_samples}},store_samples::Bool)
    # TODO what about the old store_samples_0 key
    sampler.store_samples = store_samples
end

# user_type
function set_key(sampler,::Type{Val{:user_type}},user_type::T where T)
    sampler.user_type = user_type
end

# maximum_index_set
function set_key(sampler,::Type{Val{:maximum_index_set}},maximum_index_set::Set{Index})
    # TODO only usefull in case we are running in AD mode
    sampler.maximum_index_set = maximum_index_set
end

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
