## ::::::::::::: TODO ::::::::::::::
############## 1. splitting
############## 2. do_regression (avoid expensive warm_up_samples, works hand in hand with 3.)
############## 3. repeat sampling until variance is smaller than 1.1*θ*ϵ^2 (avoids adding extra level unnecessary)
############## 4. samples0 ?
############## 5. max_index !!!
############## 6. in continuate, should only use set of keys that is currently in use!
##############    FIX: append current_index_set to Estimator
##############    let keys(estimator) return only those keys
############## 7. implement better bias formula
############## 8. PROBLEM ::: eigenfunctions are different !!!!!!!!
##############    check this with original implementation ???
##############    ===> SeparableCovarianceFunction!!!!!!!!
##############    NO! THIS IS NO SOLUTION BECAUSE OF SMOOTHNESS!!!
############## 9. add possibility for own parallel_sample! function into estimator
############## 10. MLMC with multiple
##############     make pts to interploate to a bit smaller, say 100 x 100
##############     too big! should be approx. 20 x 20!
############## 11. add rates
############## 12. catch empty bias / variance at start of loop
############## 13. merge with MonteCarloEstimator
# 14. complete history and add plotting AND reports commands
# 15. define some default cost models
#     and make ml_cost use diff
############## 16. print methods for estimator












# TODO define some standard cost models: \gamma and d, and provide multi_index_cost()
# TODO default value is Nstar = convert(Int64,ceil(32/nshifts(numberGenerators[Index(zeros(Int64,d))])))
# TODO avoid the "safety" keyword, only use doubling algorithm in QMC setting
# TODO avoid the continuate keyword, see algorithm implementation AND use the nb_of_tol keyword
# TODO how about "real" continuation ?? specify k0 and k1 (trust parameters)
# TODO what about the old store_samples_0 key
# TODO make generator states obsolete by using a new number generator at each index
# TODO QMC ::::: will just be a Matrix of nb_of_qoi x nb_of_shifts !!!!

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# 23/02: 
# then do applications: SDE, SPDE
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

#__precompile()__
module MultilevelEstimators

# load other modules
using QMC
using SpecialFunctions
using JLD
using PyPlot

# import statements
import Base: show, setindex!, run, haskey, mean, var, push!, getindex, keys, diff, γ

import QMC: ndims, nshifts, RandWrapper

# export statements
export Level, Index, diff # from index.jl

export IndexSet, SL, ML, FT, TD, HC, AD, get_index_set # from index_set.jl

export UniformMCGenerator, NormalMCGenerator, TruncatedNormalMCGenerator, UniformQMCGenerator, NormalQMCGenerator, TruncatedNormalQMCGenerator, get_point # from number_generator.jl

#export Sampler # from Sampler

#export create_estimator, MonteCarloEstimator, QuasiMonteCarloEstimator, MultiLevelMonteCarloEstimator, MultiLevelQuasiMonteCarloEstimator, MultiIndexMonteCarloEstimator, MultiIndexQuasiMonteCarloEstimator # from estimator.jl

export create_estimator#, MonteCarloEstimator, QuasiMonteCarloEstimator, MultiLevelMonteCarloEstimator, MultiLevelQuasiMonteCarloEstimator, MultiIndexMonteCarloEstimator, MultiIndexQuasiMonteCarloEstimator

export show

export run

export geometric_cost_model

# export plotting methods
export plot_E, plot_dE, plot_V, plot_dV, plot_W, plot_samples, plot_time, plot_cost

export report

export analyse

# include statements
include("index.jl")

include("index_set.jl")

include("number_generator.jl")

include("parse.jl")

include("estimator.jl")

include("history.jl")

include("sample.jl")

include("run.jl")

include("cost_models.jl")

include("monte_carlo.jl")

include("multilevel_monte_carlo.jl")

include("quasi_monte_carlo.jl")

include("multilevel_quasi_monte_carlo.jl")

include("multiindex_monte_carlo.jl")

include("print.jl")

include("plot.jl")

include("report.jl")

include("analyse.jl")

end # module
