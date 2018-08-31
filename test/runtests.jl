# runtests.jl : run all MultilevelEstimators tests

push!(LOAD_PATH,joinpath(Pkg.dir("MultilevelEstimators"),"applications","SPDE"))

using SPDE, MultilevelEstimators, QMC, Suppressor, Base.Test, JLD


using GaussianRandomFields

# basic tests
include("test_index.jl")
include("test_index_set.jl")
include("test_number_generator.jl")
include("test_estimator.jl")

## SPDE test
include("test_SPDE.jl")
include("test_analyse_SPDE.jl")
#include("generate_reports.jl")
