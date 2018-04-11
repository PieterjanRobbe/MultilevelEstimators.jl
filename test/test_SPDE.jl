## test_SPDE.jl

#=
## Monte Carlo, single qoi
@testset "Monte Carlo, single qoi      " begin
    @suppress begin
        estimator = init_lognormal_diffusion_mc()
        run(estimator,0.005)
    end
end

## Monte Carlo, multiple qoi
@testset "Monte Carlo, multiple qoi    " begin
    @suppress begin
        estimator = init_lognormal_diffusion_mc_multiple()
        run(estimator,0.005)
    end
end
=#

## Multilevel Monte Carlo, single qoi
#@testset "MLMC, single qoi             " begin
#    @suppress begin
#@show        estimator = init_lognormal_diffusion_mlmc()
#run(estimator,0.001)
#    end
#end

## Multilevel Monte Carlo, multiple qoi
#@testset "MLMC, single qoi             " begin
#    @suppress begin
@show        estimator = init_lognormal_diffusion_mlmc_multiple()
run(estimator,0.001)
#    end
#end
