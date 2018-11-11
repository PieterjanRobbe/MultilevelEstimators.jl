# runtests.jl : run all MultilevelEstimators tests

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
