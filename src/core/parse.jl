## parse.jl : parse all inputs

# parse symbol
check(settings::Dict{Symbol,T} where{T},value,symbol::Symbol) = check(settings,value,Val{symbol})

# fallback
check(::Dict{Symbol,T} where {T},::V where {V},::Type{Val{S}}) where {S} = throw(ArgumentError("could not parse key $(S)"))

# method
check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:method}}) where {V} = isa(val,AbstractIndexSet) || throw(ArgumentError("method must be of type AbstractIndexSet, got $(V)"))

# number_generator
check(settings::Dict{Symbol,T} where{T}, val::V,::Type{Val{:number_generator}}) where {V} = isa(val,NumberGenerator) || throw(ArgumentError("number_generator must be of type NumberGenerator, got $(V)"))

# sample_function
check(settings::Dict{Symbol,T} where{T}, val::V,::Type{Val{:sample_function}}) where {V} = isa(val,Function) || throw(ArgumentError("sample_function must be of type Function, got $(V)"))

# nb_of_warm_up_samples
function check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:nb_of_warm_up_samples}}) where {V}
    isa(val,Integer) || throw(ArgumentError("nb_of_warm_up_samples must be of type Integer, got $(V)"))
    ( isa(settings[:number_generator],MonteCarloNumberGenerator) && val < 2 ) && throw(ArgumentError("nb_of_warm_up_samples must be larger than 1, got $(val)"))
    ( isa(settings[:number_generator],QuasiMonteCarloNumberGenerator) && val < 1 ) && throw(ArgumentError("nb_of_warm_up_samples must be larger than 0, got $(val)"))
end

# nb_of_qoi
function check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:nb_of_qoi}}) where {V}
    isa(val,Integer) || throw(ArgumentError("nb_of_qoi must be of type Integer, got $(V)"))
    val > 0 || throw(ArgumentError("nb_of_qoi must be larger than 0, got $(val)"))
end

# continuate
function check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:continuate}}) where {V}
    isa(val,Bool) || throw(ArgumentError("continuate must be of type Bool, got $(V)"))
    if settings[:method] isa MG
        if !val
            @warn "to use a Multigrid-type estimator, set 'continuate=true'. Proceeding with 'continuate=true'..."
            settings[:continuate] = true
        end
    end
end

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
function check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:folder}}) where {V}
    isa(val,String) || throw(ArgumentError("folder must be of type String, got $(V)"))
    isdir(val) || mkpath(val)
    isempty(readdir(val)) || @warn "directory $(val) is not empty, results might be overwritten!"
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

# max_search_space
function check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:max_search_space}}) where {V}
    isa(val,AbstractIndexSet) || throw(ArgumentError("max_search_space must be of type AbstractIndexSet, got $(V)"))
    ndims(val) == ndims(settings[:method]) || throw(ArgumentError("max_search_space must be an index set in $(ndims(settings[:method])) dimensions, got $(ndims(val))"))
end

# do_regression
check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:do_regression}}) where {V} = isa(val,Bool) || throw(ArgumentError("do_regression must be of type Bool, got $(V)"))

# do_splitting
check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:do_splitting}}) where {V} = isa(val,Bool) || throw(ArgumentError("do_splitting must be of type Bool, got $(V)"))

# parallel_sample_function
check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:parallel_sample_function}}) where {V} = isa(val,Function) || throw(ArgumentError("parallel_sample_function must be of type Function, got $(V)"))

# name
check(settings::Dict{Symbol,T} where {T},val::V,::Type{Val{:name}}) where {V} = isa(val,String) || throw(ArgumentError("name must be of type String, got $(V)"))

# p0
function check(settings::Dict{Symbol,T} where{T},val::V,::Type{Val{:sample_multiplication_factor}}) where {V}
    isa(val,Real) || throw(ArgumentError("sample_multiplication_factor must be of type Real, got $(V)"))
    val >= 0 || throw(ArgumentError("p0 must be larger than or equal to 0, got $(val)"))
end
