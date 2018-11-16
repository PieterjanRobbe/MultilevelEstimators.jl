## LognormalDiffusionProblems.jl : module file
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018
module LognormalDiffusionProblems

# dependencies
using MultilevelEstimators, GaussianRandomFields

# import statements
import GaussianRandomFields: GaussianRandomFieldGenerator

# export statements
export init_lognormal

# include statements
include("init.jl")

end # module
