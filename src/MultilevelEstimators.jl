module MultilevelEstimators

# load other modules
using QMC

using FastGaussQuadrature

# import statements
import Base: +, -, *, /, ==, !=, .==

import Base: getindex, setindex!

import Base: maximum, indmax, sum, prod, one, zero

import Base.show

import Base: copy, hash, isequal

import Base.isvalid

import Base: hcat, vcat, append!

# export statements
export Index, IndexSet, createIndexSet # from Indexsetsjl

export Sampler, Settings, createSampler # from Sampler.jl

export KLexpansion, createKLexpansion, matern, compose # from GaussianFieldSamplers.jl

export MCgenerator, GaussianMCgenerator, UniformMCgenerator,
	QMCgenerator, GaussianQMCgenerator, UniformQMCgenerator, reset! # from PointSetGenerators.jl

export simulate # from MultilevelAlgorithm.jl

# include statements
include("Indexsets.jl")

include("GaussianFieldSamplers.jl")

include("PointSetGenerators.jl")

include("Sampler.jl")

include("MultilevelAlgorithm.jl")

# export all index sets
for k in instances(Kind) 
	@eval export $(symbol(k))
end 

end # module
