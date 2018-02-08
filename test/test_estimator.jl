# test_estimator.jl : test different estimator types

@testset "Estimator                    " begin

# Monte Carlo
estimator = create_estimator(method=SL(),number_generator=UniformMCGenerator(100),sample_function=(level,ξ,sampler)->geometric_brownian_motion(Level(6),ξ,sampler),tol=1e-3)
@test isa(estimator,MonteCarloEstimator)

# Quasi-Monte Carlo
estimator = create_estimator(method=SL(),number_generator=UniformQMCGenerator(100,16),sample_function=(level,ξ,sampler)->geometric_brownian_motion(Level(6),ξ,sampler),tol=1e-3)
@test isa(estimator,QuasiMonteCarloEstimator)

# Multilevel Monte Carlo
estimator = create_estimator(method=ML(),number_generator=UniformMCGenerator(100),sample_function=geometric_brownian_motion,tol=1e-3)
@test isa(estimator,MultiLevelMonteCarloEstimator)

# Multilevel Quasi-Monte Carlo
estimator = create_estimator(method=ML(),number_generator=UniformQMCGenerator(100,16),sample_function=geometric_brownian_motion,tol=1e-3)
@test isa(estimator,MultiLevelQuasiMonteCarloEstimator)

# Multi-Index Monte Carlo
estimator = create_estimator(method=TD(2),number_generator=UniformMCGenerator(100),sample_function=lognormal_diffusion,tol=1e-3)
@test isa(estimator,MultiIndexMonteCarloEstimator)

# Multi-Index Quasi-Monte Carlo
estimator = create_estimator(method=TD(2),number_generator=UniformQMCGenerator(100,16),sample_function=lognormal_diffusion,tol=1e-3)
@test isa(estimator,MultiIndexQuasiMonteCarloEstimator)

end
