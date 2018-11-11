# estimator.jl : a multilevel estimator

# TODO at some point this should be rewritten with a "BaseEstimator"
# and dedicated types for each method: avoids empty keys
struct Estimator{I<:AbstractIndexSet,G<:NumberGenerator,T<:AbstractFloat,N<:Integer,S,U,P,Q,C,H,J,D}

    # required keys
    method::I
    number_generator::G
    sample_function::Function

    # algorithm details
    nb_of_warm_up_samples::N # number of inital samples for variance estimate
    nb_of_qoi::N # number of quantities of interest
    nb_of_shifts::N # number of shifts

    # continuation
    continuate::Bool # do continuaton
    ntols::N # number of continuation steps
    p0::T # continuation parameter

    # saving options
    name::String # problem name
    folder::String # name of folder 
    store_samples::Bool # option to save samples for making pdf etc.

    # low_mem option
    # TODO low_mem::Bool # low memory option does not require storage of samples

    # internals
    samples::S # samples will be used to derive mean and variance
    samples0::S # non-difference samples for consistency check
    nsamples::P # total number of samples taken in each index
    total_work::Q # total runtime or work per index
    current_index_set::C # indices currently in use
    old_index_set::C # "old" index set in adaptive algorithm
    spill_index_set::C # "spill" indices that exceed the maximum allowed level
    max_active_set::C # maximum "active" index set to use for bias estimate
    adaptive_index_set::D # a vector with a Dict that contains all indices in the set, and their type (for plotting)
    number_generators::H # stores shifted number generator at each index

    # user_data
    has_user_data::Bool # does the estimator have user data?
    user_data::U # user_data (nothing if has_user_data is false)

    # verbose
    verbose::Bool

    # cost model
    use_cost_model::Bool
    cost_model::Function

    # conservative bias estimation (i.e., use all levels to fit rate α)
    conservative_bias_estimate::Bool

    # maximum level
    max_level::N
    max_search_space::J

    # regression instead of warm-up samples
    do_regression::Bool

    # do MSE splitting
    do_splitting::Bool

    # parallel_sample_function
    parallel_sample_function::Function

    # sample multiplication factor for QMC algorithm
    sample_multiplication_factor::T
end

const _MC = MonteCarloNumberGenerator
const _QMC = QuasiMonteCarloNumberGenerator


#const MonteCarloEstimator{T,N} = Estimator{I,G,T,N} where {I<:SL, G<:MonteCarloNumberGenerator,T,N}
#const QuasiMonteCarloEstimator{T,N} = Estimator{I,G,T,N} where {I<:SL, G<:QuasiMonteCarloNumberGenerator,T,N}
#const MultiLevelMonteCarloEstimator{T,N} = Estimator{I,G,T,N} where {I<:ML, G<:MonteCarloNumberGenerator,T,N}
#const MultiLevelQuasiMonteCarloEstimator{T,N} = Estimator{I,G,T,N} where {I<:ML, G<:QuasiMonteCarloNumberGenerator,T,N}
#const MultiIndexMonteCarloEstimator{T,N} = Estimator{I,G,T,N} where {I<:Union{TD,FT,HC,AD}, G<:MonteCarloNumberGenerator,T,N}
#const MultiIndexQuasiMonteCarloEstimator{T,N} = Estimator{I,G,T,N} where {I<:Union{TD,FT,HC,AD}, G<:QuasiMonteCarloNumberGenerator,T,N}
#const MultiGridMultiLevelMonteCarloEstimator{T,N} = Estimator{I,G,T,N} where {I<:MG{d,J} where {d,J<:ML}, G<:MonteCarloNumberGenerator,T,N}
#const MultiGridMultiLevelQuasiMonteCarloEstimator{T,N} = Estimator{I,G,T,N} where {I<:MG{d,J} where {d,J<:ML}, G<:QuasiMonteCarloNumberGenerator,T,N}
#const MultipleSemiCoarsenedMultiGridMultiIndexMonteCarloEstimator{T,N} = Estimator{I,G,T,N} where {I<:MG{d,J} where {d,J<:Union{TD,FT,HC,AD}}, G<:MonteCarloNumberGenerator,T,N}

# type aliases
#const MonteCarloTypeEstimator = Estimator{I,G,T,N} where {I,G<:MonteCarloNumberGenerator,T,N}
#const QuasiMonteCarloTypeEstimator = Estimator{I,G,T,N} where {I,G<:QuasiMonteCarloNumberGenerator,T,N}
#const SingleLevelTypeEstimator = Estimator{I} where {I<:SL}
#const MultiLevelTypeEstimator = Estimator{I} where {I<:Union{ML,MG{d,J} where {d,J<:ML}}}
#const MultiIndexTypeEstimator = Estimator{I} where {I<:Union{TD,FT,HC,AD,MG{d,J} where {d,J<:Union{TD,FT,HC,AD}}}}
#const AdaptiveMultiIndexTypeEstimator = Estimator{I} where {I<:Union{AD,MG{d,<:AD} where {d}}}
#const MultiGridMonteCarloTypeEstimator = Estimator{I,G} where {I<:MG,G<:MonteCarloNumberGenerator}
#const MultiGridQuasiMonteCarloTypeEstimator = Estimator{I,G} where {I<:MG,G<:QuasiMonteCarloNumberGenerator}

function create_estimator(;kwargs...)

    # user provided settings
    settings = Dict(kwargs)

    # check required keys
    for key in [:method, :number_generator, :sample_function]
        haskey(settings,key) || throw(ArgumentError("required key $(key) not provided"))
        check(settings,settings[key],key) # check type of arguments
    end

    # default settings
    defaults = get_default_settings(settings[:method],settings[:number_generator])

    # join the two dicts
    for (key,val) in defaults
        if !haskey(settings,key)
            settings[key] = val
        end
    end

    # check option clashes
    for (key,val) in settings
        check(settings,val,key)
    end

    # parametric types
    I = typeof(settings[:method])
    G = typeof(settings[:number_generator])
    T = Float64
    N = Int64
    S_eltype = Dict{Index{ndims(settings[:method])},Vector{T}}
    S = Matrix{S_eltype}
    U = typeof(settings[:user_data])
    P = Dict{Index{ndims(settings[:method])},N}
    Q = Dict{Index{ndims(settings[:method])},T}
    C = Set{Index{ndims(settings[:method])}}
    H = Dict{Index{ndims(settings[:method])},typeof(random_shift(settings[:number_generator]))}
    J = typeof(settings[:max_search_space])
    D = Vector{Dict{Index{ndims(settings[:method])},N}}

    # estimator internals
    m = settings[:nb_of_qoi]
    n = nb_of_shifts(settings[:number_generator])
    samples = S(undef, m, n)
    samples0 = S(undef, m, n)
    for j = 1:n
        for i = 1:m
            samples[i,j] = S_eltype()
            samples0[i,j] = S_eltype()
        end
    end
    settings[:samples] = samples
    settings[:samples0] = samples0
    settings[:nsamples] = P()
    settings[:total_work] = Q()
    settings[:has_user_data] = isa(settings[:user_data],Nothing) ? false : true
    settings[:use_cost_model] = isa(settings[:cost_model](zeros(N,ndims(settings[:method]))),Nothing) ? false : true
    settings[:current_index_set] = C()
    settings[:old_index_set] = C()
    settings[:spill_index_set] = C()
    settings[:max_active_set] = C()
    settings[:nb_of_shifts] = n
    settings[:number_generators] = H()
    settings[:adaptive_index_set] = D()

    # create estimator
    return Estimator{I,G,T,N,S,U,P,Q,C,H,J,D}([settings[name] for name in fieldnames(Estimator)]...)
end

get_default_settings(method, number_generator) = Dict{Symbol,Any}(
                                                                  :nb_of_warm_up_samples => isa(number_generator,MonteCarloNumberGenerator) ? 20 : 1,
:nb_of_qoi => 1,
:continuate => false,
:ntols => 10,
:p0 => 1.5,
:folder => "./data/",
:user_data => nothing,
:verbose => false,
:cost_model => i->nothing,
:store_samples => false,
:conservative_bias_estimate => false,
:max_level => 100,
:do_regression => true,
:do_splitting => true,
:parallel_sample_function => parallel_sample!,
:name => "",
:sample_multiplication_factor => 2,
:max_search_space => method isa AD || method isa MG{d,<:AD} where {d} ? TD(ndims(method)) : method
)

# convenience functions
haskey(estimator::Estimator,index::Index) = in(index,estimator.current_index_set)

keys(estimator::Estimator) = sort(collect(estimator.current_index_set))

#active_set(estimator::AdaptiveMultiIndexTypeEstimator) = setdiff(estimator.current_index_set,estimator.old_index_set) 

push!(estimator::Estimator,index::Index) = push!(estimator.current_index_set,index)

clear(estimator::Estimator) = begin
    if !(estimator.method isa MG)
        for index in keys(estimator)
            delete!(estimator.current_index_set,index)
            delete!(estimator.old_index_set,index)
            delete!(estimator.spill_index_set,index)
        end
        empty!(estimator.adaptive_index_set)
    end
end

ndims(estimator::Estimator) = ndims(estimator.method)

# show methods
print_name(estimator::Estimator{<:SL,<:_MC}) = "Monte Carlo estimator"
print_name(estimator::Estimator{<:ML,<:_MC}) = "Multilevel Monte Carlo estimator"
#print_name(estimator::QuasiMonteCarloEstimator) = "Quasi-Monte Carlo estimator"
#print_name(estimator::MultiLevelQuasiMonteCarloEstimator) = "Multilevel Quasi-Monte Carlo estimator"
#print_name(estimator::MultiIndexMonteCarloEstimator) = "Multi-Index Monte Carlo estimator ($(estimator.method) index set)"
#print_name(estimator::MultiIndexQuasiMonteCarloEstimator) = "Multi-Index Quasi-Monte Carlo estimator ($(estimator.method) index set)"
#print_name(estimator::MultiGridMultiLevelMonteCarloEstimator) = "MG-Multilevel Monte Carlo estimator"
#print_name(estimator::MultiGridMultiLevelQuasiMonteCarloEstimator) = "MG-Multilevel Quasi-Monte Carlo estimator"
#print_name(estimator::MultipleSemiCoarsenedMultiGridMultiIndexMonteCarloEstimator) = "MSG-Multi-Index Monte Carlo estimator ($(estimator.method) index set)"

show(io::IO, estimator::Estimator) = print(io, print_name(estimator))
