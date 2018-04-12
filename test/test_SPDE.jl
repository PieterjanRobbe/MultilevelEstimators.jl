## test_SPDE.jl

ϵ = 0.01

## Monte Carlo, single qoi
@testset "Monte Carlo, single qoi      " begin
    @suppress begin
        estimator = init_lognormal_diffusion_mc()
        run(estimator,ϵ)
    end
end

## Monte Carlo, multiple qoi
@testset "Monte Carlo, multiple qoi    " begin
    @suppress begin
        estimator = init_lognormal_diffusion_mc_multiple()
        run(estimator,ϵ)
    end
end

## Multilevel Monte Carlo, single qoi
@testset "MLMC, single qoi             " begin
    @suppress begin
        estimator = init_lognormal_diffusion_mlmc()
        run(estimator,ϵ)
    end
end

# Multilevel Monte Carlo, multiple qoi
@testset "MLMC, multiple qoi           " begin
    @suppress begin
        estimator = init_lognormal_diffusion_mlmc_multiple()
        run(estimator,ϵ)
    end
end
