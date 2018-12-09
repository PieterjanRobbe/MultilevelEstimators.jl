## internals.jl : stores estimator internals
#
# A type that stores estimator internals for specific Estimators.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## EstimatorInternals ##
abstract type AbstractEstimatorInternals end

struct EstimatorInternals{S, N, T, C, Z, T1, T2, T3} <: AbstractEstimatorInternals
    samples_diff::S
    samples::S
    nb_of_samples::N
    total_work::T
    total_time::T
    current_index_set::C
    index_set_size::Z

    sample_method_internals::T1
    index_set_internals::T2
    multigrid_internals::T3
end

function EstimatorInternals(index_set::AbstractIndexSet, sample_method::AbstractSampleMethod, settings::Dict{Symbol, Any})
    N = Int64
    T = Float64
    type_i = Index{ndims(index_set), N}
    type_d = Dict{type_i, Vector{T}}
    type_s = Matrix{type_d}
    m = settings[:nb_of_qoi]
    n = get_sample_ncols(index_set, sample_method, settings)
    type_n = Dict{type_i, N}
    type_t = Dict{type_i, T}
    type_c = Set{type_i}
    type_z = IndexSetSize{N}

    sample_method_internals = SampleMethodInternals(type_i, index_set, sample_method, settings) 
    T1 = typeof(sample_method_internals)
    index_set_internals = IndexSetInternals(type_c, type_n, index_set, sample_method, settings) 
    T2  = typeof(index_set_internals)
    multigrid_internals = MultigridInternals(index_set, T, m, n, type_t)
    T3 = typeof(multigrid_internals)

    EstimatorInternals{type_s, type_n, type_t, type_c, type_z, T1, T2, T3}(type_s(undef, m, n), type_s(undef, m, n), type_n(), type_t(), type_t(), type_c(), IndexSetSize(0, 0), sample_method_internals, index_set_internals, multigrid_internals)
end

## IndexSetSize ##
mutable struct IndexSetSize{N}
    sz::N
    max_sz::N
end

## SampleMethodInternals ##
abstract type AbstractSampleMethodInternals end

struct EmptySampleMethodInternals <: AbstractSampleMethodInternals end

struct QMCInternals{G} <: AbstractSampleMethodInternals
    generators::G
end

SampleMethodInternals(::DataType, ::AbstractIndexSet, ::AbstractSampleMethod, ::Dict{Symbol, Any}) = EmptySampleMethodInternals()

function SampleMethodInternals(type_i::DataType, index_set::AbstractIndexSet, sample_method::QMC, settings::Dict{Symbol, Any})
    type_l = typeof(ShiftedLatticeRule(settings[:point_generator]))
    type_v = Vector{type_l}
    type_g = Dict{type_i, type_v}
    QMCInternals{type_g}(type_g())
end

## IndexSetInternals ##
abstract type AbstractIndexSetInternals end

struct EmptyIndexSetInternals <: AbstractIndexSetInternals end

struct ADInternals{C, A, M} <: AbstractIndexSetInternals
    old_index_set::C
    spill_index_set::C
    boundary::C
    adaptive_index_set_log::A
    max_search_space::M
end

IndexSetInternals(::DataType, ::DataType, ::AbstractIndexSet, ::AbstractSampleMethod, ::Dict{Symbol, Any}) = EmptyIndexSetInternals()

function IndexSetInternals(type_c::DataType, type_n::DataType, index_set::AbstractAD, sample_method::AbstractSampleMethod, settings::Dict{Symbol, Any})
    type_a = Vector{type_n}
    type_m = typeof(settings[:max_search_space])
    ADInternals{type_c, type_a, type_m}(type_c(), type_c(), type_c(), type_a(undef, 0), settings[:max_search_space])
end

## MultigridInternals ##
abstract type AbstractMultigridInternals end

struct EmptyMultigridInternals <: AbstractMultigridInternals end

struct MultigridInternals{A, W} <: AbstractMultigridInternals
    accumulator::A
    weights::W
end

MultigridInternals(::AbstractIndexSet, ::DataType, ::Integer, ::Integer, ::DataType) = EmptyMultigridInternals()

function MultigridInternals(::MG, T::DataType, m::Integer, n::Integer, type_t::DataType)
    type_a = Matrix{Vector{T}}

    MultigridInternals{type_a, type_t}(type_a(undef, m, n), type_t())
end

# dispatch on the size of the Matrix that holds the samples
get_sample_ncols(index_set::AbstractIndexSet, ::MC, settings::Dict{Symbol, Any}) = 1

function get_sample_ncols(index_set::AbstractIndexSet, ::QMC, settings::Dict{Symbol, Any})
    nbshifts = settings[:nb_of_shifts]
    L = settings[:max_index_set_param]
    idx_set = get_max_index_set(index_set, settings, L)
    maximum(nbshifts.(idx_set))
end

get_max_index_set(index_set::AbstractIndexSet, settings::Dict{Symbol, Any}, L::Integer) = get_index_set(index_set, L)
get_max_index_set(index_set::AbstractAD, settings::Dict{Symbol, Any}, L::Integer) = get_index_set(settings[:max_search_space], L)
