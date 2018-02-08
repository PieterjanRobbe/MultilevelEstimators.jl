# runtests.jl : run all MultilevelEstimators tests

push!(LOAD_PATH,"../applications/SDE")
push!(LOAD_PATH,"../applications/SPDE")

using QMC
using MultilevelEstimators
using SDE
using SPDE
using Base.Test

include("test_index.jl")
include("test_index_set.jl")
include("test_number_generator.jl")
include("test_estimator.jl")
