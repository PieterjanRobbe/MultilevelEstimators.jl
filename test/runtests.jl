# runtests.jl : run all MultilevelEstimators tests
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for
# Multilevel Monte Carlo Methods

using Distributed, MultilevelEstimators, Printf, Suppressor, Test

## unit tests ##

include("index.jl")

include("index_set.jl")

include("distribution.jl")

include("lattice.jl")

include("estimator.jl")

include("regression.jl")
