## parse.jl : parse all inputs and set the corresponding sampler options

# parse symbol
setindex!(sampler,value,symbol::Symbol) = setindex!(sampler,value,Val{symbol})

# fallback option
setindex!(::Sampler, ::Any, ::Type{Val{T}}) where T = throw(ArgumentError("could not parse key $(key)!"))

## REQUIRED ##

# method
setindex!(sampler,method::IndexSet,::Type{Val{:method}}) = nothing

# number_generator
setindex!(sampler,number_generator::NumberGenerator,::Type{Val{:number_generator}}) = nothing

# sample_functiom
function setindex!(sampler,sample_function::Function,::Type{Val{:sample_function}})
	sampler.sample_function = sample_function
end

# ml_sample_function
function setindex!(sampler,ml_sample_function,::Type{Val{:ml_sample_function}})
	sampler.sample_function = ml_sample_function
end

# user_data
setindex!(sampler,user_data,::Type{Val{:user_data}}) = nothing

## OPTIONAL ##

# cost_model
function setindex!(sampler,cost_model::Function,::Type{Val{:cost_model}})
    sampler.use_sample_run_time == true && warn("conflicting options use_sample_run_time and cost_model detected, using cost_model")
    sampler.use_sample_run_time = false
    sampler.cost_model = cost_model
end

# use_sample_run_time
function setindex!(sampler,use_sample_run_time::Bool,::Type{Val{:use_sample_run_time}})
	( use_sample_run_time && !( typeof(sampler.cost_model) <: EmptyFunction ) ) && warn("conflicting options use_sample_run_time and cost_model detected, using use_sample_run_time")
    sampler.use_sample_run_time = use_sample_run_time
	sampler.cost_model = use_sample_run_time ? EmptyFunction() : sampler.cost_model
end

# nb_of_qoi
function setindex!(sampler,nb_of_qoi::N where {N<:Integer},::Type{Val{:nb_of_qoi}})
    nb_of_qoi > 0 || throw(ArgumentError("nb_of_qoi must be positive, got $(nb_of_qoi)"))
    sampler.nb_of_qoi = nb_of_qoi
end

# nb_of_warm_up_samples
function setindex!(sampler,nb_of_warm_up_samples::N where {N<:Integer},::Type{Val{:nb_of_warm_up_samples}})
    nb_of_warm_up_samples > 0 || throw(ArgumentError("nb_of_warm_up_samples must be positive, got $(nb_of_warm_up_samples)!"))
    sampler.nb_of_warm_up_samples = nb_of_warm_up_samples
end

# nb_of_tol
function setindex!(sampler,nb_of_tol::N where {N<:Integer},::Type{Val{:nb_of_tol}})
    nb_of_tol > 0 || throw(ArgumentError("nb_of_tol must be positive, got $(nb_of_tol)!"))
    sampler.nb_of_tol = nb_of_tol
end

# maximum_level
function setindex!(sampler,maximum_level::Level,::Type{Val{:maximum_level}})
    typeof(sampler.method) <: ML || throw(ArgumentError("cannot set maximum level if sampler uses a $(sampler.method)"))
    isempty(sampler.valid_indices) || warn("conflicting options maximum_level and valid_indices detected, using maximum_level")
    sampler.valid_indices = get_indexset(ML(),maximum_level[1])
end

# maximum_index
function setindex!(sampler,maximum_index::Index{d},::Type{Val{:maximum_level}}) where {d}
    typeof(sampler.method) <: Union{SL,ML} && throw(ArgumentError("cannot set maximum index if sampler uses a $(sampler.method)"))
    ndims(sampler.method) == d || throw(DimensionMismatch("number of dimensions in maximum_index does not agree with number of dimensions of the sampler (expected $(ndims(sampler.method)), got $(d))"))
    isempty(sampler.valid_indices) || warn("conflicting options maximum_index and valid_indices detected, using maximum_index")
    N = maximum(maximum_index)
    idxset = FT(ndims(sampler.method),δ=1./maximum_index*N)
    sampler.valid_indices = get_indexset(idxset,N)
end

# valid_indices
function setindex!(sampler,::Type{Val{:valid_indices}},valid_indices::Set{Index{d}}) where {d}
    ndims(sampler.method) == d || throw(DimensionMismatch("number of dimensions in valid_indices does not agree with number of dimensions of the sampler (expected $(ndims(sampler.method)), got $(d))"))
    isempty(sampler.valid_indices) || warn("setting valid_indices, but already specified")
    is_valid_index_set(valid_indices) || throw(ArgumentError("valid_indices is not a valid index set"))
    sampler.valid_indices = valid_indices
end

# low_mem
function setindex!(sampler,low_mem::Bool,::Type{Val{:low_mem}})
    sampler.low_mem = low_mem
end
