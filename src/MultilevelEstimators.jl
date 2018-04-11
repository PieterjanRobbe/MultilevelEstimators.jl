# TODO define some standard cost models: \gamma and d, and provide multi_index_cost()
# TODO define use_sample_run_time = true (default), unless a cost_model is provided
# TODO default value is Nstar = convert(Int64,ceil(32/nshifts(numberGenerators[Index(zeros(Int64,d))])))
# TODO avoid the "safety" keyword, only use doubling algorithm in QMC setting
# TODO avoid the continuate keyword, see algorithm implementation AND use the nb_of_tol keyword
# TODO how about "real" continuation ?? specify k0 and k1 (trust parameters)
# TODO what about the old store_samples_0 key
# TODO make generator states obsolete by using a new number generator at each index
# TODO make dicts as functions: E(sampler,index), V(sampler, index); E(sampler), V(sampler)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# 23/02: 
# estimator input problem: immutable struct
# then do applications: SDE, SPDE
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

#__precompile()__
module MultilevelEstimators

# load other modules
using QMC
using SpecialFunctions
using JLD
using PyPlot

# import statements
import Base: show, setindex!, run, haskey, mean, var, push!, getindex, keys, diff, γ

import QMC: ndims, nshifts, RandWrapper

# export statements
export Level, Index, diff # from index.jl

export IndexSet, SL, ML, FT, TD, HC, AD, get_index_set # from index_set.jl

export UniformMCGenerator, NormalMCGenerator, TruncatedNormalMCGenerator, UniformQMCGenerator, NormalQMCGenerator, TruncatedNormalQMCGenerator, get_point # from number_generator.jl

#export Sampler # from Sampler

#export create_estimator, MonteCarloEstimator, QuasiMonteCarloEstimator, MultiLevelMonteCarloEstimator, MultiLevelQuasiMonteCarloEstimator, MultiIndexMonteCarloEstimator, MultiIndexQuasiMonteCarloEstimator # from estimator.jl

export create_estimator#, MonteCarloEstimator, QuasiMonteCarloEstimator, MultiLevelMonteCarloEstimator, MultiLevelQuasiMonteCarloEstimator, MultiIndexMonteCarloEstimator, MultiIndexQuasiMonteCarloEstimator

export show

export run

# include statements
include("index.jl")

include("index_set.jl")

include("number_generator.jl")

#include("sampler.jl")

include("parse.jl")

include("estimator.jl")

include("history.jl")

include("sample.jl")

include("run.jl")

include("monte_carlo.jl")

include("multilevel_monte_carlo.jl")

end # module
