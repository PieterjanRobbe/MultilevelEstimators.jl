## LognormalDiffusionProblems.jl : module file
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018
module LognormalDiffusionProblems

# dependencies
using Distributed, FFTW, GaussianRandomFields, Interpolations, MultilevelEstimators, NotSoSimpleMultigrid, PaddedViews, SimpleMultigrid, Statistics

# import statements
import GaussianRandomFields: GaussianRandomFieldGenerator

import SimpleMultigrid: MultigridIterable

# export statements
export init_lognormal, sample_lognormal, Qoi1, Qoi2, Qoi3, MGSolver, MSGSolver

# include statements
include("init.jl")

include("sample.jl")

include("FMG.jl")

end # module
