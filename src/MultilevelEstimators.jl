module MultilevelEstimators

# load other modules
using QMC

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

import QMC: getPoint

# export statements
export Index, SL, ML, FT, TD, HC, AD, getIndexSet # from Indexsets.jl

export Sampler, setup, save, load, reset # from Sampler.jl

export KLExpansion, MaternKernel, compose # from GaussianFieldSamplers.jl

export GaussianMCgenerator, UniformMCgenerator, GaussianQMCgenerator, UniformQMCgenerator, reset # from PointSetGenerators.jl

export simulate # from MultilevelAlgorithm.jl

# include statements
include("Indexsets.jl")

include("GaussianFieldSamplers.jl")

include("PointSetGenerators.jl")

include("Sampler.jl")

include("MultilevelAlgorithm.jl")

end # module
