## test_SPDE.jl

ϵ₁= 0.1
ϵ₂= 0.005

## Monte Carlo, single qoi
@testset "MC, single qoi               " begin
    @suppress begin
        estimator = init_SPDE_mc(continuate=true)
        run(estimator,ϵ₁)
    end
end

## Monte Carlo, multiple qoi
@testset "MC, multiple qoi             " begin
    @suppress begin
        estimator = init_SPDE_mc_multiple(continuate=true)
        run(estimator,ϵ₁)
    end
end

## Multilevel Monte Carlo, single qoi
@testset "MLMC, single qoi             " begin
    @suppress begin
        estimator = init_SPDE_mlmc(continuate=true)
        run(estimator,ϵ₂)
    end
end

# Multilevel Monte Carlo, multiple qoi
@testset "MLMC, multiple qoi           " begin
    @suppress begin
        estimator = init_SPDE_mlmc_multiple(continuate=true)
        run(estimator,ϵ₂)
    end
end

#=
## Quasi-Monte Carlo, single qoi
@testset "QMC, single qoi              " begin
    @suppress begin
        estimator = init_SPDE_qmc(continuate=true)
        run(estimator,ϵ₁)
    end
end

# Quasi-Monte Carlo, multiple qoi
@testset "QMC, multiple qoi            " begin
    @suppress begin
        estimator = init_SPDE_qmc_multiple(continuate=true)
        run(estimator,ϵ₁)
    end
end

## Multilevel Quasi-Monte Carlo, single qoi
@testset "MLQMC, single qoi            " begin
    @suppress begin
        estimator = init_SPDE_mlqmc(continuate=true)
        run(estimator,ϵ₂)
    end
end

# Multilevel Quasi-Monte Carlo, multiple qoi
@testset "MLQMC, multiple qoi          " begin
    @suppress begin
        estimator = init_SPDE_mlqmc_multiple(continuate=true)
        run(estimator,ϵ₂)
    end
end

# Multi-Index Monte Carlo, single qoi
@testset "MIMC, single qoi             " begin
    @suppress begin
        estimator = init_SPDE_mimc(continuate=true)
        run(estimator,ϵ₂)
    end
end

# Multi-Index Monte Carlo, multiple qoi
@testset "MIMC, multiple qoi           " begin
    @suppress begin
        estimator = init_SPDE_mimc_multiple(continuate=true)
        run(estimator,ϵ₂)
    end
end

# Multi-Index Quasi-Monte Carlo, single qoi
@testset "MIQMC, single qoi            " begin
    @suppress begin
        estimator = init_SPDE_miqmc(continuate=true)
        run(estimator,ϵ₂)
    end
end

# Multi-Index Quasi-Monte Carlo, multiple qoi
@testset "MIQMC, multiple qoi          " begin
    @suppress begin
        estimator = init_SPDE_miqmc_multiple(continuate=true)
        run(estimator,ϵ₂)
    end
end

# Adaptive Multi-Index Monte Carlo, single qoi
@testset "AMIMC, single qoi            " begin
    @suppress begin
        estimator = init_SPDE_amimc(continuate=true)
        run(estimator,ϵ₂)
    end
end

# Adaptive Multi-Index Monte Carlo, multiple qoi
@testset "AMIMC, multiple qoi          " begin
    @suppress begin
        estimator = init_SPDE_amimc_multiple(continuate=true)
        run(estimator,ϵ₂)
    end
end

# Adaptive Multi-Index Quasi-Monte Carlo, single qoi
@testset "AMIQMC, single qoi           " begin
    @suppress begin
        estimator = init_SPDE_amiqmc(continuate=true)
        run(estimator,ϵ₂)
    end
end

# Adaptive Multi-Index Quasi-Monte Carlo, multiple qoi
@testset "AMIQMC, multiple qoi         " begin
    @suppress begin
        estimator = init_SPDE_amiqmc_multiple(continuate=true)
        run(estimator,ϵ₂)
    end
end

## Multigrid Multilevel Monte Carlo, single qoi
@testset "MG-MLMC, single qoi          " begin
    @suppress begin
        estimator = init_SPDE_mgmlmc(continuate=true)
        run(estimator,ϵ₂)
    end
end

# Multigrid Multilevel Monte Carlo, multiple qoi
@testset "MG-MLMC, multiple qoi        " begin
    @suppress begin
        estimator = init_SPDE_mgmlmc_multiple(continuate=true)
        run(estimator,ϵ₂)
    end
end

## Multiple Semi-Coarsened Multigrid Multi-Index Monte Carlo, single qoi
@testset "MSG-MIMC, single qoi         " begin
    @suppress begin
        estimator = init_SPDE_msgmimc(continuate=true)
        run(estimator,ϵ₂)
    end
end

## Multiple Semi-Coarsened Multigrid Multi-Index Monte Carlo, multiple qoi
@testset "MSG-MIMC, multiple qoi       " begin
    @suppress begin
        estimator = init_SPDE_msgmimc_multiple(continuate=true)
        run(estimator,ϵ₂)
    end
end

## Multiple Semi-Coarsened Multigrid Adaptive Multi-Index Monte Carlo, single qoi
@testset "MSG-MIMC, single qoi         " begin
    @suppress begin
        estimator = init_SPDE_msgamimc(continuate=true)
        run(estimator,ϵ₂)
    end
end

## Multiple Semi-Coarsened Multigrid Adaptive Multi-Index Monte Carlo, multiple qoi
@testset "MSG-MIMC, multiple qoi       " begin
    @suppress begin
        estimator = init_SPDE_msgamimc_multiple(continuate=true)
        run(estimator,ϵ₂)
    end
end

## Multigrid Multilevel Quasi-Monte Carlo, single qoi
@testset "MG-MLQMC, single qoi         " begin
    @suppress begin
        estimator = init_SPDE_mgmlqmc(continuate=true)
        run(estimator,ϵ₂)
    end
end

# Multigrid Multilevel Quasi-Monte Carlo, multiple qoi
@testset "MG-MLQMC, multiple qoi       " begin
    @suppress begin
        estimator = init_SPDE_mgmlqmc_multiple(continuate=true)
        run(estimator,ϵ₂)
    end
end
=#
