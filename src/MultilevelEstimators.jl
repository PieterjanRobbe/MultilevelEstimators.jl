## MultilevelEstimators.jl : module file
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

module MultilevelEstimators

## dependencies ##

using Colors, Dates, DelimitedFiles, Distributed, JLD2, LinearAlgebra, Printf, Random, SpecialFunctions, Statistics

## import statements ##

import Base: show, diff, getindex, filter, ndims, eltype, run, push!, keys, zero, isless, >, <, ≥, ≤, +, length, one

import Random: GLOBAL_RNG

import Statistics: mean, var

## export statements ##

export Level, Index

export AbstractIndexSet, SL, ML, FT, TD, HC, AD, MG, get_index_set

export AbstractDistribution, Uniform, Normal, TruncatedNormal, Weibull 

export AbstractSampleMethod, MC, QMC

export Estimator

export LatticeRule32

export run

export report

## include statements ##

# core

include("core/check_inputs.jl")

include("core/index.jl")

include("core/index_sets.jl")

include("core/distributions.jl")

include("core/sample_methods.jl")

include("core/options.jl")

include("core/internals.jl")

include("core/estimators.jl")

include("core/parsers.jl")

include("core/lattices.jl")

include("core/print.jl")

include("core/sample.jl")

include("core/history.jl")

include("core/run.jl")

# methods

include("methods/monte_carlo.jl")

include("methods/multilevel_monte_carlo.jl")

include("methods/multiindex_monte_carlo.jl")

# utils

include("utils/report.jl")

include("utils/tex.jl")

end # module
