## LognormalDiffusionProblems.jl : module file
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018
module LognormalDiffusionProblems

# dependencies
using Distributed, FFTW, GaussianRandomFields, Interpolations, MultilevelEstimators, NotSoSimpleMultigrid, PaddedViews, Random, Retry, SimpleMultigrid, SpecialFunctions, Statistics

# import statements
import GaussianRandomFields: GaussianRandomFieldGenerator

import SimpleMultigrid: MultigridIterable, MultigridCycle

# export statements
export init_lognormal, sample_lognormal, Qoi1, Qoi2, Qoi3, Qoi4, MGSolver, MSGSolver, AnalyzeV, AnalyzeFMG, NonIsotropicMatern, compute_grf

# include statements
include("init.jl")

include("sample.jl")

include("FMG.jl")

include("NonIsotropicMatern.jl")

end # module
