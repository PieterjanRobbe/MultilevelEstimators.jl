# runtests.jl : run all MultilevelEstimators tests

push!(LOAD_PATH,joinpath(Pkg.dir("MultilevelEstimators"),"applications","SPDE"))

using QMC, MultilevelEstimators, SPDE, Suppressor, Base.Test, JLD

estimator = init_SPDE_mlmc(continuate=true)

#@show SPDE.SPDE_sample_mg_single(Level(5), randn(100_000), estimator.user_data)
#@show SPDE.SPDE_sample_mg_single(Level(5), randn(100_000), estimator.user_data)
#@show SPDE.SPDE_sample_mg_single(Level(5), randn(100_000), estimator.user_data)
#@show SPDE.SPDE_sample_mg_single(Level(5), randn(100_000), estimator.user_data)
@show SPDE.SPDE_sample_single(Level(5), randn(100_000), estimator.user_data)
SPDE.SPDE_sample_multiple(Level(5), randn(100_000), estimator.user_data)

#=
wp = CachingPool(workers())
f(i) = estimator.sample_function(Level(5),randn(100_000),estimator.user_data)
t = @elapsed all_samples = pmap(wp,f,1:10)
@show typeof(all_samples)
=#

#=
# basic tests
include("test_index.jl")
include("test_index_set.jl")
include("test_number_generator.jl")
include("test_estimator.jl")

## SPDE test
include("test_SPDE.jl")
include("test_analyse_SPDE.jl")
#include("generate_reports.jl")
=#
