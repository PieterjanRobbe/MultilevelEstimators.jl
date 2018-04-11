## sampler.jl : defintion of a Sampler

struct Sampler{I<:IndexSet,G<:NumberGenerator,U,V,D1,D2}

    ## REQUIRED ##
    index_set::I
    number_generator::G
    data::U	

    sample_function::Function
    ml_sample_function::Function

    ## OPTIONAL ##
    cost_model::Function
    use_sample_run_time::Bool
    nb_of_qoi::Int64
    nb_of_warm_up_samples::Int64
    valid_indices::V
    low_mem::Bool

    ## INTERNAL ##
    samples::D1
    diff_samples::D1
    run_times::D2

end

## nullables
struct EmptyFunction <: Function end
struct EmptyData end

## basic sampler
function Sampler(method::IndexSet{d},number_generator::G,user_data::U) where {d,G<:NumberGenerator,U}
    Sampler{typeof(method),G,U,Vector{Index{d}},Dict{Index{d},Array{Float64,3}},Dict{Index{d},Float64}}(
        method, # index_set
        number_generator, # number_generator
        user_data, # data
        EmptyFunction(), # sample_function
        ml_sample, # ml_sample_function
        EmptyFunction(), # cost_model
        false, # use_sample_run_time
        1, # nb_of_qoi
        20, # nb_of_warm_up_samples
        Index{d}[], # valid_indices
        false, # low_mem
        Dict{Index{d},Array{Float64,3}}(), # samples
        Dict{Index{d},Array{Float64,3}}(), # diff_samples
        Dict{Index{d},Float64}() # run_time
    )
end

## utilities

#--------------- sampling -----------------

function ml_sample(index::Index{d} where {d},sampler::Sampler)
    return false
end

function sample(index::Index{d} where {d},nsamples::N where {N<:Integer},sampler::Sampler)
    return false
end

#--------------- inspection -----------------

function MSE(sampler::Sampler)
    return false
end

function E(sampler::Sampler,index::Index{d} where {d})
    return false
end

function V(sampler::Sampler,index::Index{d} where {d})
    return false
end

function W(sampler::Sampler,index::Index{d} where {d})
    return false
end
