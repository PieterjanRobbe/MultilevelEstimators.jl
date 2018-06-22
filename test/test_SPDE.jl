## test_SPDE.jl

ϵ₁= 0.1
ϵ₂= 0.005

## Monte Carlo, single qoi
@testset "MC, single qoi               " begin
    @suppress begin
        estimator = init_lognormal_diffusion_mc()
        run(estimator,ϵ₁)
    end
end

## Monte Carlo, multiple qoi
@testset "MC, multiple qoi             " begin
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

# Multi-Index Monte Carlo, single qoi
@testset "MIMC, single qoi             " begin
    @suppress begin
        estimator = init_lognormal_diffusion_mimc()
        run(estimator,ϵ₂)
    end
end

# Multi-Index Monte Carlo, multiple qoi
@testset "MIMC, multiple qoi           " begin
    @suppress begin
        estimator = init_lognormal_diffusion_mimc_multiple()
        run(estimator,ϵ₂)
    end
end

# Multi-Index Quasi-Monte Carlo, single qoi
@testset "MIQMC, single qoi            " begin
    @suppress begin
        estimator = init_lognormal_diffusion_miqmc()
        run(estimator,ϵ₂)
    end
end

# Multi-Index Quasi-Monte Carlo, multiple qoi
@testset "MIQMC, multiple qoi          " begin
    @suppress begin
        estimator = init_lognormal_diffusion_miqmc_multiple()
        run(estimator,ϵ₂)
    end
end
