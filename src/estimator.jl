# estimator.jl : a multilevel estimator

## Estimator ##
struct Estimator{S,H,T}
    sampler::S
    history::H
    tol::T
end

"""
    create_estimator(kwargs...)

Main method to create a multilevel estimator.

# Examples
```jldoctest
julia>

```
"""
function create_estimator(;kwargs...)
    settings = Dict(kwargs)

    # required
    for key in [:method, :number_generator, :tol]
        haskey(settings,key) || throw(ArgumentError("required key $(key) not provided"))
    end

    # create basic sampler
    data = haskey(settings,:user_data) ? settings[:user_data] : EmptyData()
    sampler = Sampler(settings[:method],settings[:number_generator],data) 

    # read tolerances
    if ( typeof(settings[:tol]) <: Vector{T} where {T<:Real} || typeof(settings[:tol]) <: T where {T<:Real} )
        tol = settings[:tol]
    else
        throw(ArgumentError("incorrect type of required parameter tol"))
    end

    # set sample function
    if haskey(settings,:sample_function)
        sampler[:sample_function] = settings[:sample_function]
    elseif haskey(settings,:ml_sample_function)
        sampler[:ml_sample_function] = settings[:ml_sample_function]
    else
        throw(ArgumentError("required key sample_function not provided")) 
    end

    # remove all required keys keys
    for req_key in [:method, :number_generator, :user_data, :tol]
        delete!(settings,req_key)
    end
        
    # set optional keys
    for (key, value) in settings
        setindex!(sampler, value, key)
    end

    # cost model or sample_run_time
    if !haskey(settings,:cost_model) && 
        ( !haskey(settings,:use_sample_run_time) || !settings[:use_sample_run_time]) 
        setindex!(sampler, true, :use_sample_run_time)
    end

    # history
    history = History()

    # estimator
    Estimator(sampler,history,tol)
end

## Type aliases ##

##### HERE #####
## TODO
const MonteCarloEstimator = Estimator{Sampler{SL,NumberGenerator{D,MC}}} where {D}
const QuasiMonteCarloEstimator = Estimator{Sampler{SL,NumberGenerator{D,QMC}}} where {D}
const MultiLevelMonteCarloEstimator = Estimator{Sampler{I,NumberGenerator{D,G}}} where {I<:ML,D,G<:MC}
MultilevelEstimators.Estimator{S,H,T} where {S<:MultilevelEstimators.Sampler{I,MultilevelEstimators.NumberGenerator{D,G}} where {I<:ML,D,G<:MultilevelEstimators.MC},H,T}
const MultiLevelQuasiMonteCarloEstimator = Estimator{Sampler{ML,NumberGenerator{D,QMC}}} where {D}
const MultiIndexMonteCarloEstimator = Estimator{Sampler{I,NumberGenerator{D,MC}}} where {D,I<:Union{FT,TD,HC,AD}}
const MultiIndexQuasiMonteCarloEstimator = Estimator{Sampler{I,NumberGenerator{D,QMC}}} where {D,I<:Union{FT,TD,HC,AD}}

