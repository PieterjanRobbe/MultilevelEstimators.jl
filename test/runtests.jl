## runtests.jl : run all test files

using MultilevelEstimators
using Base.Test

# toggle printing mode
verbose = true

# test indexsets
include("test_index_sets.jl")
include("test_number_generators.jl")
#include("test_multilevelalgorithm.jl")
