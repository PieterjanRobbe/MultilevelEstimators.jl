# test_estimator.jl : test different estimator types

@testset "Estimator                    " begin

# Monte Carlo
@suppress estimator = create_estimator(
    method=SL(),
    number_generator=UniformMCGenerator(100),
    sample_function=i->nothing
)
@test isa(estimator,MultilevelEstimators.MonteCarloEstimator)

# Quasi-Monte Carlo
@suppress estimator = create_estimator(
    method=SL(),
    number_generator=UniformQMCGenerator(100,16),
    sample_function=i->nothing
)
@test isa(estimator,MultilevelEstimators.QuasiMonteCarloEstimator)

# Multilevel Monte Carlo
@suppress estimator = create_estimator(
    method=ML(),
    number_generator=UniformMCGenerator(100),
    sample_function=i->nothing
)
@test isa(estimator,MultilevelEstimators.MultiLevelMonteCarloEstimator)

# Multilevel Quasi-Monte Carlo
@suppress estimator = create_estimator(
    method=ML(),
    number_generator=UniformQMCGenerator(100,20),
    sample_function=i->nothing
)
@test isa(estimator,MultilevelEstimators.MultiLevelQuasiMonteCarloEstimator)

# Multi-Index Monte Carlo
@suppress estimator = create_estimator(
    method=TD(2),
    number_generator=UniformMCGenerator(100),
    sample_function=i->nothing
)
@test isa(estimator,MultilevelEstimators.MultiIndexMonteCarloEstimator)

# Multi-Index Quasi-Monte Carlo
@suppress estimator = create_estimator(
    method=TD(2),
    number_generator=UniformQMCGenerator(100,8),
    sample_function=i->nothing
)
@test isa(estimator,MultilevelEstimators.MultiIndexQuasiMonteCarloEstimator)

# test other options
@suppress estimator = create_estimator(
    method=SL(),
    number_generator=UniformMCGenerator(100),
    sample_function=i->nothing,
    verbose=true,
    continuate=true,
    folder="some_folder",
    nb_of_qoi=120,
    cost_model=i->prod(2.^i)
)
@test isa(estimator,MultilevelEstimators.MonteCarloEstimator)

@test_throws ArgumentError create_estimator(
    method=TD(2),
    number_generator=UniformQMCGenerator(100,16),
)

@test_throws ArgumentError create_estimator(
    number_generator=UniformQMCGenerator(100,165),
    sample_function=i->nothing,
)

@test_throws ArgumentError create_estimator(
    method=TD(2),
    sample_function=i->nothing,
)

@suppress @test_throws ArgumentError create_estimator(
    method=TD(2),
    number_generator=UniformQMCGenerator(100,23),
    sample_function=i->nothing,
    continuate=1.
)

@test_throws ArgumentError create_estimator(
    method=TD(2),
    number_generator=UniformMCGenerator(100),
    sample_function=i->nothing,
    folder=1
)

# etc

end
