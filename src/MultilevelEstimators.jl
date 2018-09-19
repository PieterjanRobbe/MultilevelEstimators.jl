__precompile__()
module MultilevelEstimators

# load other modules
using QMC, SpecialFunctions, JLD#, PyPlot

# import statements
import Base: show, setindex!, run, haskey, mean, var, push!, getindex, keys, diff, γ

import QMC: ndims, nshifts, RandWrapper

# export statements
export Level, Index, diff # from index.jl

export IndexSet, SL, ML, FT, TD, HC, AD, MG, get_index_set # from index_set.jl

export UniformMCGenerator, NormalMCGenerator, TruncatedNormalMCGenerator, UniformQMCGenerator, NormalQMCGenerator, TruncatedNormalQMCGenerator, get_point # from number_generator.jl

export create_estimator, show, clear # from estimator.jl

export run # from run.jl

export geometric_cost_model # from cost_models.jl

export save # from history.jl

export report # from report.jl

export plot_E, plot_dE, plot_V, plot_dV, plot_W, plot_samples, plot_time, plot_cost, plot_time_vs_rmse # from plot.jl

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

include("methods/multigrid_multilevel_monte_carlo.jl")

include("methods/multiple_semicoarsened_multigrid_multiindex_monte_carlo.jl")

include("utils/plot.jl")

include("utils/tex.jl")

include("utils/report.jl")

end # module
