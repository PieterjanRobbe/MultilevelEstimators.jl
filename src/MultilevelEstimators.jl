## MultilevelEstimators.jl : module file
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

module MultilevelEstimators

## dependencies ##

using Colors, Dates, DelimitedFiles, Distributed, JLD2, LinearAlgebra, Printf, Random, SpecialFunctions, Statistics

## import statements ##

import Base: filter, ndims, eltype, show, diff, sum, :, print, length, getindex, run, push!, isless, >, <, ≥, ≤, +, zero, one, keys 

import Random: GLOBAL_RNG

import Statistics: mean, var

## export statements ##

export Level, Index

export AbstractIndexSet, SL, ML, FT, TD, HC, ZC, AD, U, MG, get_index_set

export AbstractDistribution, Uniform, Normal, TruncatedNormal, Weibull , transform

export AbstractSampleMethod, MC, QMC

export LatticeRule32, get_point

export Estimator

#export run

#export report

## include statements ##

# core

include("core/check_input.jl")

include("core/index.jl")

include("core/index_set.jl")

include("core/distribution.jl")

include("core/sample_method.jl")

include("core/lattice.jl")

include("core/estimator.jl")

include("core/options.jl")

include("core/parse.jl")

include("core/internals.jl")

#include("core/print.jl")

#include("core/sample.jl")

#include("core/history.jl")

#include("core/run.jl")

# methods

#include("methods/monte_carlo.jl")

#include("methods/multilevel_monte_carlo.jl")

#include("methods/multiindex_monte_carlo.jl")

# utils

#include("utils/report.jl")

#include("utils/tex.jl")

end # module
