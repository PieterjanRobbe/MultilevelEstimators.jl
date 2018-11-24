## LognormalDiffusionProblems.jl : module file
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018
module LognormalDiffusionProblems

# dependencies
using MultilevelEstimators, GaussianRandomFields, SimpleMultigrid, PaddedViews, Interpolations

# import statements
import GaussianRandomFields: GaussianRandomFieldGenerator

import SimpleMultigrid: MultigridIterable, μ_cycle!, norm_of_residu, P

# export statements
export init_lognormal, sample_lognormal, Qoi1, Qoi2, Qoi3

# include statements
include("init.jl")

include("sample.jl")

include("FMG.jl")

end # module
