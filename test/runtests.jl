using MultilevelEstimators
using Base.Test

# toggle printing mode
verbose = true

# test indexsets
include("test_indexsets.jl")
include("test_multilevelalgorithm.jl")
