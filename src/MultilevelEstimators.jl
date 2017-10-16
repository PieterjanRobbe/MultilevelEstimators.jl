module MultilevelEstimators

# load other modules
using   ArgParse,
        FastGaussQuadrature, 
        Interpolations, 
        NLsolve, 
        PmapProgressMeter,
        ProgressMeter, 
        Reexport,
        Retry, 
        SpecialFunctions 

@reexport using QMC

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

#import Base: ndims

import QMC: ndims, getPoint, nshifts

# export statements
export Index, ndims # from indices.jl

export IndexSet, SL, ML, FT, TD, HC, AD, get_index_set, get_boundary, pretty_print, is_admissable # from index_set.jl

export GaussianMCgenerator, UniformMCgenerator, GaussianQMCgenerator, UniformQMCgenerator, reset # from number_generators.jl

export CovarianceFunction, MaternCovarianceFunction # from random_fields.jl

export EmptyGaussianRandomFieldSampler # from random_field_samplers.jl

### ### ###

export Sampler, setup, save, load, reset, sample, ndims, mi_cost  # from Sampler.jl

export KLExpansion, MaternKernel, compose, preload_eigenfunctions # from GaussianFieldSamplers.jl


export simulate # from MultilevelAlgorithm.jl
export prashant_simulate # from prashant_algorithm.jl

export FEsolve, sample, sort # FEsolve only for test.jl

export @debug

export analyse

# TRM
export leastSquaresFit
#

# include statements
include("indices.jl")

include("index_sets.jl")

include("number_generators.jl")

include("random_fields.jl")

include("random_field_samplers.jl")

### ### ###

#include("karhunen_loeve_expansion.jl")

include("sampler.jl")

include("parse_options.jl")

include("parse_sampler.jl")

include("multigrid_rng.jl")

include("multigrid_sampler.jl")
include("giles_multigrid_sampler.jl")
include("prashant_algorithm.jl")

#include("MultilevelAlgorithm.jl")

#include("analyse.jl")


end # module
