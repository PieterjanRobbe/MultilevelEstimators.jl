# runtests.jl : run all MultilevelEstimators tests

using Test, Printf, Random, Distributed
using MultilevelEstimators

## unit tests ##
include("index.jl")
include("index_sets.jl")
include("distributions.jl")
include("estimators.jl")
include("lattices.jl")

## lognormal tests ##
push!(LOAD_PATH,joinpath(dirname(pathof(MultilevelEstimators)), "..", "applications", "lognormal"))
using LognormalDiffusionProblems, GaussianRandomFields
include("lognormal.jl")
