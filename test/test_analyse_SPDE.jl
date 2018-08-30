## test_analyse_SPDE.jl

## Analyse SPDE, multilevel case
@testset "analyse SPDE, ML case        " begin
    @suppress begin
        estimator = init_SPDE_analyse_ml()
        analyse(estimator, nsamples=100)
    end
end

#=
## Analyse SPDE, multi-index case
@testset "analyse SPDE, MI case        " begin
    @suppress begin
        estimator = init_SPDE_analyse_mi()
        analyse(estimator, nsamples=100)
    end
end
=#
