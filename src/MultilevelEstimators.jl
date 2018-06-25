# TODO add possibility for own parallel_sample! function into estimator
# TODO define some standard cost models: \gamma and d, and provide multi_index_cost()
# TODO default value is Nstar = convert(Int64,ceil(32/nshifts(numberGenerators[Index(zeros(Int64,d))])))
# TODO how about "real" continuation ?? specify k0 and k1 (trust parameters)
# TODO what about the old store_samples_0 key

#__precompile()__
module MultilevelEstimators

# load other modules
using QMC, SpecialFunctions, JLD, PyPlot

# import statements
import Base: show, setindex!, run, haskey, mean, var, push!, getindex, keys, diff, γ

import QMC: ndims, nshifts, RandWrapper

# export statements
export Level, Index, diff # from index.jl

export IndexSet, SL, ML, FT, TD, HC, AD, get_index_set # from index_set.jl

export UniformMCGenerator, NormalMCGenerator, TruncatedNormalMCGenerator, UniformQMCGenerator, NormalQMCGenerator, TruncatedNormalQMCGenerator, get_point # from number_generator.jl

export create_estimator, show # from estimator.jl

export run # from run.jl

export geometric_cost_model # from cost_models.jl

# export plotting methods
export save # from history.jl

export report # from report.jl

export plot_E, plot_dE, plot_V, plot_dV, plot_W, plot_samples, plot_time, plot_cost # from plot.jl

export analyse # from analyse.jl

# include statements
include("core/index.jl")

include("core/index_set.jl")

include("core/number_generator.jl")

include("core/parse.jl")

include("core/estimator.jl")

include("core/history.jl")

include("core/sample.jl")

include("core/run.jl")

include("core/cost_models.jl")

include("core/print.jl")

include("core/analyse.jl")

include("methods/monte_carlo.jl")

include("methods/multilevel_monte_carlo.jl")

include("methods/quasi_monte_carlo.jl")

include("methods/multilevel_quasi_monte_carlo.jl")

include("methods/multiindex_monte_carlo.jl")

include("methods/multiindex_quasi_monte_carlo.jl")

include("methods/adaptive_multiindex_monte_carlo.jl")

include("utils/plot.jl")

include("utils/tex.jl")

include("utils/report.jl")

end # module
