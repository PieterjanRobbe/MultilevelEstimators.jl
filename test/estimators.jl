## estimators.jl : unit testing for estimators.jl
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

@testset "Estimator                    " begin
    index_set = ML()
    ad_index_set = AD(3)
    mc_sample_method = MC()
    qmc_sample_method = QMC()
    sample_function = () -> ()
    distribution = Normal()
    distributions = [Normal(), Normal(), Normal()]
    type_mlmc = Estimator{typeof(index_set), typeof(mc_sample_method)}
    type_mlqmc = Estimator{typeof(index_set), typeof(qmc_sample_method)}
    type_amlqmc = Estimator{typeof(ad_index_set), typeof(qmc_sample_method)}

    @test Estimator(index_set, mc_sample_method, sample_function, distributions) isa type_mlmc
    @test Estimator(index_set, qmc_sample_method, sample_function, distribution) isa type_mlqmc

    # some valid estimators here ...
    
    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; nb_of_warm_up_samples=false)
    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; nb_of_warm_up_samples=0)
    @test Estimator(index_set, mc_sample_method, sample_function, distribution; nb_of_warm_up_samples=10) isa type_mlmc
    
    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; nb_of_qoi=false)
    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; nb_of_qoi=0)
    @test Estimator(index_set, mc_sample_method, sample_function, distribution; nb_of_qoi=10) isa type_mlmc
    
    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; continuate="blabla")
    @test Estimator(index_set, mc_sample_method, sample_function, distribution; continuate=false) isa type_mlmc
    
    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; nb_of_tols="blabla")
    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; nb_of_tols=0)
    @test Estimator(index_set, mc_sample_method, sample_function, distribution; nb_of_tols=20) isa type_mlmc
    
    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; continuation_mul_factor="blabla")
    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; continuation_mul_factor=NaN)
    @test Estimator(index_set, qmc_sample_method, sample_function, distribution; continuation_mul_factor=2) isa type_mlqmc

    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; folder=1.0)
    @test Estimator(index_set, mc_sample_method, sample_function, distribution; folder=".") isa type_mlmc

    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; name=1.0)
    @test Estimator(index_set, mc_sample_method, sample_function, distribution; name="MyEstimator") isa type_mlmc

    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; save_samples=1)
    @test Estimator(index_set, mc_sample_method, sample_function, distribution; save_samples=true) isa type_mlmc

    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; verbose=1)
    @test Estimator(index_set, mc_sample_method, sample_function, distribution; verbose=true) isa type_mlmc

    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; cost_model="blabla")
    @test Estimator(index_set, mc_sample_method, sample_function, distribution; cost_model=i->2^sum(i)) isa type_mlmc

    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; robustify_bias_estimate="blabla")
    @test Estimator(index_set, mc_sample_method, sample_function, distribution; robustify_bias_estimate=false) isa type_mlmc

    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; do_mse_splitting=1)
    @test Estimator(index_set, mc_sample_method, sample_function, distribution; do_mse_splitting=false) isa type_mlmc

    @test Estimator(index_set, mc_sample_method, sample_function, distribution; min_splitting=0.6) isa type_mlmc
    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; min_splitting=2)

    @test Estimator(index_set, mc_sample_method, sample_function, distribution; max_splitting=0.8) isa type_mlmc
    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; max_splitting=0.4)
    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; max_splitting=0.6, min_splitting=0.7)

    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; max_index_set_param=1.0)
    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; max_index_set_param=-1)
    @test Estimator(index_set, mc_sample_method, sample_function, distribution; max_index_set_param=6) isa type_mlmc
    
    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; sample_mul_factor="blabla")
    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; sample_mul_factor=Inf)
    @test Estimator(index_set, qmc_sample_method, sample_function, distribution; sample_mul_factor=1.5) isa type_mlqmc
    @test Estimator(index_set, qmc_sample_method, sample_function, distribution; sample_mul_factor=2) isa type_mlqmc

    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; nb_of_workers=Inf)
    @test Estimator(index_set, mc_sample_method, sample_function, distribution; nb_of_workers=nworkers()) isa type_mlmc
    @test Estimator(index_set, mc_sample_method, sample_function, distribution; nb_of_workers=i->2) isa type_mlmc

    @test_throws ArgumentError Estimator(index_set, qmc_sample_method, sample_function, distribution; nb_of_shifts="blabla")
    @test Estimator(index_set, qmc_sample_method, sample_function, distribution; nb_of_shifts=16) isa type_mlqmc
    @test Estimator(index_set, qmc_sample_method, sample_function, distribution; nb_of_shifts=i->32) isa type_mlqmc 

    @test_throws ArgumentError Estimator(index_set, qmc_sample_method, sample_function, distribution; point_generator=3)

    @test_throws ArgumentError Estimator(index_set, mc_sample_method, sample_function, distribution; do_regression="blabla")

    @test_throws ArgumentError Estimator(AD(2), mc_sample_method, sample_function, distribution; max_search_space=true)
    @test_throws ArgumentError Estimator(ad_index_set, mc_sample_method, sample_function, distribution; max_search_space=TD(2))
    @test Estimator(ad_index_set, qmc_sample_method, sample_function, distribution; max_search_space=FT(3)) isa type_amlqmc
end
