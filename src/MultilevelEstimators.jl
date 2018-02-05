#__precompile()__
module MultilevelEstimators

# load other modules
using QMC
using SpecialFunctions

# import statements
import Base.show

import QMC: ndims, nshifts, RandWrapper

# export statements
export Level, Index # from index.jl

export SL, ML, FT, TD, HC, AD, get_indexset # from index_set.jl

export UniformMCGenerator, NormalMCGenerator, TruncatedNormalMCGenerator, UniformQMCGenerator, NormalQMCGenerator, TruncatedNormalQMCGenerator, get_point

export show

# include statements
include("index.jl")

include("index_set.jl")

include("number_generator.jl")

end # module
