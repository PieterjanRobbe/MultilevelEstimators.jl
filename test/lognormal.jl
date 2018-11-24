## lognormal.jl : testing for LognormalDiffusionProblems.jl
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

# tmp folder where history files will be stored
my_temp_folder = tempdir()

## Monte Carlo ##
@testset "(lognormal) MC               " begin
    estimator = init_lognormal(SL(), MC(), 
                               nb_of_coarse_dofs=8,
                               grf_generator=KarhunenLoeve(100),
                               covariance_function=Matern(1,2),
                               folder=my_temp_folder
                              )
    run(estimator, 1e-2)
    estimator = init_lognormal(SL(), MC(), 
                               nb_of_coarse_dofs=16,
                               qoi=Qoi2(),
                               cost_model=(level)->16*2^level[1],
                               folder=my_temp_folder
                              )
    run(estimator, 1e-1)
    estimator = init_lognormal(SL(), MC(),
                               qoi=Qoi3(),
                               nb_of_qoi=16,
                               nb_of_coarse_dofs=32,
                               folder=my_temp_folder
                              )
    run(estimator, 1e-2)
end

## Multilevel Monte Carlo ##
@testset "(lognormal) MLMC             " begin
    estimator = init_lognormal(ML(), MC(), 
                               nb_of_coarse_dofs=8,
                               grf_generator=KarhunenLoeve(100),
                               covariance_function=Matern(1,2),
                               max_index_set_param=3,
                               folder=my_temp_folder
                              )
    run(estimator, 1e-3)
    estimator = init_lognormal(ML(), MC(), 
                               nb_of_coarse_dofs=16,
                               qoi=Qoi2(),
                               cost_model=(level)->16*2^level[1],
                               max_index_set_param=2,
                               continuate=false,
                               do_mse_splitting=false,
                               do_regression=false,
                               folder=my_temp_folder
                              )
    run(estimator, 1e-1)
    estimator = init_lognormal(ML(), MC(),
                               qoi=Qoi3(),
                               nb_of_qoi=16, 
                               nb_of_coarse_dofs=32,
                               max_index_set_param=3,
                               folder=my_temp_folder
                              )
    run(estimator, 1e-2)
end
