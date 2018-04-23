## test_SPDE.jl

ϵ₁= 0.05
ϵ₂= 0.005

## Monte Carlo, single qoi
@testset "Monte Carlo, single qoi      " begin
    @suppress begin
        estimator = init_lognormal_diffusion_mc()
        run(estimator,ϵ₁)
    end
end

## Monte Carlo, multiple qoi
@testset "Monte Carlo, multiple qoi    " begin
    @suppress begin
        estimator = init_lognormal_diffusion_mc_multiple()
        run(estimator,ϵ₁)
    end
end

## Multilevel Monte Carlo, single qoi
@testset "MLMC, single qoi             " begin
    @suppress begin
        estimator = init_lognormal_diffusion_mlmc()
        run(estimator,ϵ₂)
    end
end

# Multilevel Monte Carlo, multiple qoi
@testset "MLMC, multiple qoi           " begin
    @suppress begin
        estimator = init_lognormal_diffusion_mlmc_multiple()
        run(estimator,ϵ₂)
    end
end

## Quasi-Monte Carlo, single qoi
@testset "QMC, single qoi              " begin
    @suppress begin
        estimator = init_lognormal_diffusion_qmc()
        run(estimator,ϵ₁)
    end
end

# Quasi-Monte Carlo, multiple qoi
@testset "QMC, multiple qoi            " begin
    @suppress begin
        estimator = init_lognormal_diffusion_qmc_multiple()
        run(estimator,ϵ₁)
    end
end

## Multilevel Quasi-Monte Carlo, single qoi
@testset "MLQMC, single qoi            " begin
    @suppress begin
        estimator = init_lognormal_diffusion_mlqmc()
        run(estimator,ϵ₂)
    end
end

# Multilevel Quasi-Monte Carlo, multiple qoi
@testset "MLQMC, multiple qoi          " begin
    @suppress begin
        estimator = init_lognormal_diffusion_mlqmc_multiple()
        run(estimator,ϵ₂)
    end
end


