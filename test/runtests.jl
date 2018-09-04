# runtests.jl : run all MultilevelEstimators tests

push!(LOAD_PATH,joinpath(Pkg.dir("MultilevelEstimators"),"applications","SPDE"))

using SPDE, MultilevelEstimators, QMC, Suppressor, Base.Test, JLD

estimator = init_SPDE_msgmimc()
#estimator = init_SPDE_mgmlmc()

samples = [estimator.sample_function(Index(2,2),randn(500),estimator.user_data) for i in 1:1000]

@show mean(first.(samples))
@show mean(last.(samples))

#run(estimator,0.001)

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
