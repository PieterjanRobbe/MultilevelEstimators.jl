module MultilevelEstimators

# load other modules
using Reexport

@reexport using QMC

using Interpolations, FastGaussQuadrature, SpecialFunctions, Retry, NLsolve, ArgParse

using ProgressMeter
using PmapProgressMeter

# import statements
import Base: +, -, *, .*, /, ==, !=, .==

import Base: getindex, setindex!, length

import Base: maximum, indmax

import Base: sum, prod, mean, std, diff

import Base: one, zero

import Base: show, reset, sort

import Base: copy, hash, isequal, isless

import Base: hcat, vcat, append!

import Base: collect

import Base: ndims

import QMC: getPoint, nshifts

# export statements
export Index, SL, ML, FT, TD, HC, AD, get_index_set, get_boundary, pretty_print, is_admissable # from index_set.jl

export Sampler, setup, save, load, reset, sample, ndims, mi_cost  # from Sampler.jl

export KLExpansion, MaternKernel, compose, preload_eigenfunctions # from GaussianFieldSamplers.jl

export GaussianMCgenerator, UniformMCgenerator, GaussianQMCgenerator, UniformQMCgenerator, PseudoRNG, reset # from number_generators.jl

export simulate # from MultilevelAlgorithm.jl
export prashant_simulate # from prashant_algorithm.jl

export FEsolve, sample, sort # FEsolve only for test.jl

export @debug

export analyse

# TRM
export leastSquaresFit
#

# include statements
include("parse_options.jl")

include("index_sets.jl")

include("GaussianFieldSamplers.jl")

include("number_generators.jl")
include("multigrid_rng.jl")

include("Sampler.jl")
include("multigrid_sampler.jl")
include("giles_multigrid_sampler.jl")
include("prashant_algorithm.jl")

include("MultilevelAlgorithm.jl")

include("analyse.jl")


end # module
