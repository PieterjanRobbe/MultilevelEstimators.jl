## MultilevelEstimators.jl : module file
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for 
# Multilevel Monte Carlo Methods

module MultilevelEstimators

## dependencies ##

using Dates, DelimitedFiles, Distributed, JLD2, LinearAlgebra, Printf, Random, SpecialFunctions, Statistics

## import statements ##

import Base: filter, ndims, eltype, show, diff, sum, :, print, length, getindex, run, push!, isless, >, <, ≥, ≤, +, zero, one, keys, haskey

import Random: GLOBAL_RNG

import Statistics: mean, var

## export statements ##

export Level, Index

export AbstractIndexSet, SL, ML, FT, TD, HC, ZC, AD, U, get_index_set, is_admissable, length

export AbstractDistribution, Uniform, Normal, TruncatedNormal, Weibull , transform

export AbstractSampleMethod, MC, QMC

export LatticeRule32, ShiftedLatticeRule, get_point

export Estimator, mean, var

export History, haskey

export run

## include statements ##

# core

include("check_input.jl")

include("index.jl")

include("index_set.jl")

include("distribution.jl")

include("sample_method.jl")

include("lattice.jl")

include("estimator.jl")

include("options.jl")

include("parse.jl")

include("internals.jl")

include("print.jl")

include("sample.jl")

include("history.jl")

include("run.jl")

include("methods.jl")

end # module
