## estimator.jl : unit testing for estimator.jl
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for
# Multilevel Monte Carlo Methods

@testset "Estimator                    " begin
	ml = ML()
	ad = AD(3)
	mc = MC()
	qmc = QMC()
	qoi = () -> ()
	distr = Normal()
	distrs = [Normal() for i in 1:1000]

	@test Estimator(ml, mc, qoi, distrs) isa Estimator{<:ML, <:MC}
	@test Estimator(ml, qmc, qoi, distr) isa Estimator{<:ML, <:QMC}

    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; nb_of_warm_up_samples=false)
    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; nb_of_warm_up_samples=0)
    @test Estimator(ml, mc, qoi, distr; nb_of_warm_up_samples=10) isa Estimator{<:ML, <:MC}
    
    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; nb_of_qoi=false)
    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; nb_of_qoi=0)
    @test Estimator(ml, mc, qoi, distr; nb_of_qoi=10) isa Estimator{<:ML, <:MC}
    
    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; continuate="blabla")
    @test Estimator(ml, mc, qoi, distr; continuate=false) isa Estimator{<:ML, <:MC}
    
    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; nb_of_tols="blabla")
    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; nb_of_tols=0)
    @test Estimator(ml, mc, qoi, distr; nb_of_tols=20) isa Estimator{<:ML, <:MC}
    
    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; continuation_mul_factor="blabla")
    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; continuation_mul_factor=NaN)
	@test Estimator(ml, qmc, qoi, distr; continuation_mul_factor=2) isa Estimator{<:ML, <:QMC}

    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; folder=1.0)
    @test Estimator(ml, mc, qoi, distr; folder=".") isa Estimator{<:ML, <:MC}

    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; name=1.0)
    @test Estimator(ml, mc, qoi, distr; name="MyEstimator") isa Estimator{<:ML, <:MC}

    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; save_samples=1)
    @test Estimator(ml, mc, qoi, distr; save_samples=true) isa Estimator{<:ML, <:MC}

    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; verbose=1)
    @test Estimator(ml, mc, qoi, distr; verbose=true) isa Estimator{<:ML, <:MC}

    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; cost_model="blabla")
    @test Estimator(ml, mc, qoi, distr; cost_model=i->2^sum(i)) isa Estimator{<:ML, <:MC}

    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; robustify_bias_estimate="blabla")
    @test Estimator(ml, mc, qoi, distr; robustify_bias_estimate=false) isa Estimator{<:ML, <:MC}

    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; do_mse_splitting=1)
    @test Estimator(ml, mc, qoi, distr; do_mse_splitting=false) isa Estimator{<:ML, <:MC}

    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; min_splitting=2)
    @test Estimator(ml, mc, qoi, distr; min_splitting=0.6) isa Estimator{<:ML, <:MC}

    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; max_splitting=0.4)
    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; max_splitting=0.6, min_splitting=0.7)
    @test Estimator(ml, mc, qoi, distr; max_splitting=0.8) isa Estimator{<:ML, <:MC}

    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; max_index_set_param=1.0)
    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; max_index_set_param=-1)
    @test Estimator(ml, mc, qoi, distr; max_index_set_param=6) isa Estimator{<:ML, <:MC}
    
    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; sample_mul_factor="blabla")
    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; sample_mul_factor=Inf)
    @test Estimator(ml, qmc, qoi, distr; sample_mul_factor=1.5) isa Estimator{<:ML, <:QMC}
    @test Estimator(ml, qmc, qoi, distr; sample_mul_factor=2) isa Estimator{<:ML, <:QMC}

    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; nb_of_workers=Inf)
    @test Estimator(ml, mc, qoi, distr; nb_of_workers=nworkers()) isa Estimator{<:ML, <:MC}
    @test Estimator(ml, mc, qoi, distr; nb_of_workers=i->2) isa Estimator{<:ML, <:MC}

    @test_throws ArgumentError Estimator(ml, qmc, qoi, distr; nb_of_shifts="blabla")
    @test Estimator(ml, qmc, qoi, distr; nb_of_shifts=16) isa Estimator{<:ML, <:QMC}
    @test Estimator(ml, qmc, qoi, distr; nb_of_shifts=i->32) isa Estimator{<:ML, <:QMC}

    @test_throws ArgumentError Estimator(ml, qmc, qoi, distr; point_generator=3)

    @test_throws ArgumentError Estimator(ml, mc, qoi, distr; do_regression="blabla")

    @test_throws ArgumentError Estimator(AD(2), mc, qoi, distr; max_search_space=true)
    @test_throws ArgumentError Estimator(ad, mc, qoi, distr; max_search_space=TD(2))
    @test Estimator(ad, qmc, qoi, distr; max_search_space=FT(3)) isa Estimator{<:AD, <:QMC}

	# TODO add more tests... 
end
