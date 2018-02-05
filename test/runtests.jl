# runtests.jl : run all MultilevelEstimators tests

using QMC
using MultilevelEstimators
using Base.Test

include("test_index.jl")
include("test_index_set.jl")
include("test_number_generator.jl")
