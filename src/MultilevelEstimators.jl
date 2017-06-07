module MultilevelEstimators

# load other modules
using Reexport

@reexport using QMC

using FastGaussQuadrature

# import statements
import Base: +, -, *, .*, /, ==, !=, .==

import Base: getindex, setindex!, length

import Base: maximum, indmax

import Base: sum, prod

import Base: one, zero

import Base: show, reset, sort

import Base: copy, hash, isequal, isless

import Base: hcat, vcat, append!

import Base: collect

import Base: ndims

import QMC: getPoint, nshifts

# export statements
export Index, SL, ML, FT, TD, HC, AD, getIndexSet, getBoundary # from Indexsets.jl

export Sampler, setup, save, load, reset, sample, ndims # from Sampler.jl

export KLExpansion, MaternKernel, compose # from GaussianFieldSamplers.jl

export GaussianMCgenerator, UniformMCgenerator, GaussianQMCgenerator, UniformQMCgenerator, reset # from PointSetGenerators.jl

export simulate # from MultilevelAlgorithm.jl

export FEsolve, sample, sort # FEsolve only for test.jl

export @debug

# include statements
include("ParseArgs.jl")

include("Indexsets.jl")

include("GaussianFieldSamplers.jl")

include("PointSetGenerators.jl")

include("Sampler.jl")

include("MultilevelAlgorithm.jl")

end # module
