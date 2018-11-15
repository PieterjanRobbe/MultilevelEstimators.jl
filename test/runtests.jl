# runtests.jl : run all MultilevelEstimators tests

using Test, Printf, Random, Distributed
using MultilevelEstimators

# unit tests
include("index.jl")
include("index_sets.jl")
include("distributions.jl")
include("estimators.jl")
include("lattices.jl")

#=
import MultilevelEstimators
push!(LOAD_PATH,joinpath(dirname(pathof(MultilevelEstimators)),"..","applications","SPDE"))

using SPDE, MultilevelEstimators, QMC, Suppressor, Test, JLD, Printf

# basic tests
include("test_index.jl")
include("test_index_set.jl")
include("test_number_generator.jl")
include("test_estimator.jl")

## SPDE test
include("test_SPDE.jl")
include("test_analyse_SPDE.jl")
#include("generate_reports.jl")
=#
