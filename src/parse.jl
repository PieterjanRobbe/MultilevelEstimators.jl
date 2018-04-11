## parse.jl : parse all inputs and set the corresponding sampler options

# parse symbol
check(settings::Dict{Symbol,T} where{T},value,symbol::Symbol) = check(settings,value,Val{symbol})

# fallback
check(::Dict{Symbol,T} where {T},::V where {V},::Type{Val{S}}) where {S} = throw(ArgumentError("could not parse key $(S)"))

# method
check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:method}}) where {V} = isa(val,IndexSet) || throw(ArgumentError("method must be of type IndexSet, got $(V)"))

# number_generator
check(settings::Dict{Symbol,T} where{T}, val::V,::Type{Val{:number_generator}}) where {V} = isa(val,NumberGenerator) || throw(ArgumentError("number_generator must be of type NumberGenerator, got $(V)"))

# sample_function
check(settings::Dict{Symbol,T} where{T}, val::V,::Type{Val{:sample_function}}) where {V} = isa(val,Function) || throw(ArgumentError("sample_function must be of type Function, got $(V)"))

# nb_of_warm_up_samples
function check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:nb_of_warm_up_samples}}) where {V}
    isa(val,Integer) || throw(ArgumentError("nb_of_warm_up_samples must be of type Integer, got $(V)"))
    val > 0 || throw(ArgumentError("nb_of_warm_up_samples must be larger than 0, got $(val)"))
end

# nb_of_qoi
function check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:nb_of_qoi}}) where {V}
    isa(val,Integer) || throw(ArgumentError("nb_of_qoi must be of type Integer, got $(V)"))
    val > 0 || throw(ArgumentError("nb_of_qoi must be larger than 0, got $(val)"))
end

# continuate
check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:continuate}}) where {V} = isa(val,Bool) || throw(ArgumentError("continuate must be of type Bool, got $(V)"))

# ntols
function check(settings::Dict{Symbol,T} where{T},val::V,::Type{Val{:ntols}}) where {V}
    isa(val,Integer) || throw(ArgumentError("ntols must be of type Integer, got $(V)"))
    val > 0 || throw(ArgumentError("ntols must be larger than 0, got $(val)"))
end

# p0
function check(settings::Dict{Symbol,T} where{T},val::V,::Type{Val{:p0}}) where {V}
    isa(val,Real) || throw(ArgumentError("p0 must be of type Real, got $(V)"))
    val > 0 || throw(ArgumentError("p0 must be larger than 0, got $(val)"))
end

# folder
function check(settings::Dict{Symbol,T} where{T},val::V,::Type{Val{:folder}}) where {V}
    isa(val,String) || throw(ArgumentError("folder must be of type String, got $(V)"))
    isdir(val) || mkpath(val)
    isempty(readdir(val)) || warn("directory $(val) is not empty, results might be overwritten!")
    settings[:folder] = val[end] == '/' ? val : val*"/"
end

check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:store_samples}}) where {V} = isa(val,Bool) || throw(ArgumentError("store_samples must be of type Bool, got $(V)"))

# low_mem
# TODO check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:low_mem}}) where {V} = isa(val,Bool) || throw(ArgumentError("low_mem must be of type Bool, got $(V)"))

# user_data
check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:user_data}}) where {V} = nothing

# continuate
check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:verbose}}) where {V} = isa(val,Bool) || throw(ArgumentError("verbose must be of type Bool, got $(V)"))

# cost_model
check(settings::Dict{Symbol,T} where{T}, val::V,::Type{Val{:cost_model}}) where {V} = isa(val,Function) || throw(ArgumentError("cost_model must be of type Function, got $(V)"))

# conservative_bias_estimate
check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:conservative_bias_estimate}}) where {V} = isa(val,Bool) || throw(ArgumentError("conservative_bias_estimate must be of type Bool, got $(V)"))

# max_level
function check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:max_level}}) where {V}
    isa(val,Integer) || throw(ArgumentError("max_level must be of type Integer, got $(V)"))
    val < 0 && throw(ArgumentError("max_level must be larger than or equal to 0, got $(val)"))
end

# do_regression
check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:do_regression}}) where {V} = isa(val,Bool) || throw(ArgumentError("do_regression must be of type Bool, got $(V)"))

#=
# parse symbol
setindex!(sampler::Sampler,value,symbol::Symbol) = setindex!(sampler,value,Val{symbol})

# fallback option
setindex!(::Sampler, ::V where {V}, ::Type{Val{T}}) where {T} = throw(ArgumentError("could not parse key $(T)"))

## REQUIRED ##

# sample_functiom
function setindex!(sampler::Sampler,sample_function::Function,::Type{Val{:sample_function}})
    sampler.sample_function = sample_function
end

# ml_sample_function
function setindex!(sampler::Sampler,ml_sample_function::Function,::Type{Val{:ml_sample_function}})
    sampler.sample_function = ml_sample_function
end

## OPTIONAL ##

# cost_model
function setindex!(sampler::Sampler,cost_model::Function,::Type{Val{:cost_model}})
    sampler.use_sample_run_time == true && warn("conflicting options use_sample_run_time and cost_model detected, using cost_model")
    sampler.use_sample_run_time = false
    sampler.cost_model = cost_model
end

# use_sample_run_time
function setindex!(sampler::Sampler,use_sample_run_time::Bool,::Type{Val{:use_sample_run_time}})
    ( use_sample_run_time && !( typeof(sampler.cost_model) <: EmptyFunction ) ) && warn("conflicting options use_sample_run_time and cost_model detected, using use_sample_run_time")
    sampler.use_sample_run_time = use_sample_run_time
    sampler.cost_model = use_sample_run_time ? EmptyFunction() : sampler.cost_model
end

# nb_of_qoi
function setindex!(sampler::Sampler,nb_of_qoi::N where {N<:Integer},::Type{Val{:nb_of_qoi}})
    nb_of_qoi > 0 || throw(ArgumentError("nb_of_qoi must be positive, got $(nb_of_qoi)"))
    sampler.nb_of_qoi = nb_of_qoi
end

# nb_of_warm_up_samples
function setindex!(sampler::Sampler,nb_of_warm_up_samples::N where {N<:Integer},::Type{Val{:nb_of_warm_up_samples}})
    nb_of_warm_up_samples > 0 || throw(ArgumentError("nb_of_warm_up_samples must be positive, got $(nb_of_warm_up_samples)!"))
    sampler.nb_of_warm_up_samples = nb_of_warm_up_samples
end

# maximum_level
function setindex!(sampler::Sampler,maximum_level::Level,::Type{Val{:maximum_level}})
    typeof(sampler.index_set) <: ML || throw(ArgumentError("cannot set maximum level if sampler uses a $(sampler.index_set), use maximum_index instead"))
    maximum_level[1] > 1 || throw(ArgumentError("I need at least two levels in maximum_level, got $(maximum_level[1])"))
    isempty(sampler.valid_indices) || warn("conflicting options maximum_level and valid_indices detected, using maximum_level")
    sampler.valid_indices = get_index_set(ML(),maximum_level[1])
end

# maximum_index
function setindex!(sampler::Sampler,maximum_index::Index{d},::Type{Val{:maximum_index}}) where {d}
    typeof(sampler.index_set) <: Union{SL,ML} && throw(ArgumentError("cannot set maximum_index if sampler uses a $(sampler.index_set), use maximum_level instead"))
    ndims(sampler.index_set) == d || throw(DimensionMismatch("number of dimensions in maximum_index does not agree with number of dimensions of the sampler (expected $(ndims(sampler.index_set)), got $(d))"))
    all(maximum_index .> 1) || throw(ArgumentError("I need at least two levels in each direction in maximum_index, got $(maximum_index)"))
    isempty(sampler.valid_indices) || warn("conflicting options maximum_index and valid_indices detected, using maximum_index")
    N = maximum(maximum_index)
    idxset = FT(ndims(sampler.index_set),δ=collect(1./maximum_index.*N))
    sampler.valid_indices = get_index_set(idxset,N)
end

# valid_indices
function setindex!(sampler::Sampler,valid_indices::Vector{I} where {I<:Index{d}},::Type{Val{:valid_indices}}) where {d}
    ndims(sampler.index_set) == d || throw(DimensionMismatch("number of dimensions in valid_indices does not agree with number of dimensions of the sampler (expected $(ndims(sampler.index_set)), got $(d))"))
    isempty(sampler.valid_indices) || warn("setting valid_indices, but already specified")
    is_valid_index_set(valid_indices) || throw(ArgumentError("valid_indices is not a valid index set"))
    sampler.valid_indices = valid_indices
end

# low_mem
# TODO fix: make this in create method, pick appropriate Sampler(...) with different types for samples, samples_diff
function setindex!(sampler::Sampler,low_mem::Bool,::Type{Val{:low_mem}})
    sampler.low_mem = low_mem
end
=#
