# runtests.jl : run all MultilevelEstimators tests

push!(LOAD_PATH,Pkg.dir("MultilevelEstimators/applications/SPDE"))

using QMC, MultilevelEstimators, SPDE, Suppressor, Base.Test

# basic tests
include("test_index.jl")
include("test_index_set.jl")
include("test_number_generator.jl")
include("test_estimator.jl")

## SPDE test
include("test_SPDE.jl")
include("test_analyse_SPDE.jl")
