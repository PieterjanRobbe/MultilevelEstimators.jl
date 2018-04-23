## test_analyse_SPDE.jl

## Analyse SPDE, multilevel case
@testset "analyse SPDE, ML case        " begin
    @suppress begin
        estimator = init_lognormal_diffusion_analyse_ml()
        analyse(estimator)
    end
end
