## MultilevelEstimators.jl : module file
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018
module MultilevelEstimators

# dependencies
using LinearAlgebra, SpecialFunctions, Distributed, Random, DelimitedFiles, Printf, Dates, Statistics

# import statements
import Base: show, diff, getindex, filter, ndims, eltype, run, push!, keys, zero, isless, >, <, ≥, ≤, +

import Random: GLOBAL_RNG

import Statistics: mean, var

# export statements
export Level, Index

export AbstractIndexSet, SL, ML, FT, TD, HC, AD, MG, get_index_set

export AbstractDistribution, Uniform, Normal, TruncatedNormal, Weibull 

export AbstractSampleMethod, MC, QMC

export Estimator

export LatticeRule32

export run

# include statements
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

include("core/run.jl")

include("methods/monte_carlo.jl")

include("methods/multilevel_monte_carlo.jl")

#=

# load other modules
using Printf, Statistics, Distributed, Dates, SpecialFunctions, JLD, QMC, LinearAlgebra, SparseArrays

# import statements
import Base: show, setindex!, run, haskey, push!, getindex, keys, diff #, γ
import Statistics: mean, var

import QMC: ndims, nshifts, RandWrapper

# export statements
export Level, Index, diff # from index.jl

export AbstractIndexSet, SL, ML, FT, TD, HC, AD, MG, get_index_set # from index_set.jl

export UniformMCGenerator, NormalMCGenerator, TruncatedNormalMCGenerator, UniformQMCGenerator, NormalQMCGenerator, TruncatedNormalQMCGenerator, get_point # from number_generator.jl

export create_estimator, show, clear # from estimator.jl

export run # from run.jl

export geometric_cost_model # from cost_models.jl

export save # from history.jl

export report # from report.jl

export plot_E, plot_dE, plot_V, plot_dV, plot_W, plot_samples, plot_time, plot_cost, plot_time_vs_rmse # from plot.jl

export analyse # from analyse.jl

# include statements
include("core/index.jl")

include("core/index_set.jl")

include("core/number_generator.jl")

include("core/parse.jl")

include("core/estimator.jl")

include("core/history.jl")

include("core/sample.jl")

include("core/run.jl")

include("core/cost_models.jl")

include("core/print.jl")

include("core/analyse.jl")

include("methods/monte_carlo.jl")

include("methods/multilevel_monte_carlo.jl")

#include("methods/quasi_monte_carlo.jl")

#include("methods/multilevel_quasi_monte_carlo.jl")

#include("methods/multiindex_monte_carlo.jl")

#include("methods/multiindex_quasi_monte_carlo.jl")

#include("methods/multigrid_multilevel_monte_carlo.jl")

#include("methods/multiple_semicoarsened_multigrid_multiindex_monte_carlo.jl")

#include("methods/multigrid_multilevel_quasi_monte_carlo.jl")

#include("utils/plot.jl")

include("utils/tex.jl")

include("utils/report.jl")
=#

end # module
