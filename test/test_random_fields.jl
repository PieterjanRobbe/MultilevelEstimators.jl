## test_random_fields.jl : tests for random_fields.jl

## MaternCovarianceFunction ##
verbose && print("testing matern covariance function...")

m = MaternCovarianceFunction(2,0.1,1.,2.5,1.)
@test typeof(m) <: MaternCovarianceFunction
m = MaternCovarianceFunction(2,1,1.,2.5,1.)
@test typeof(m) <: MaternCovarianceFunction
m = MaternCovarianceFunction(2,0.1,1,2.5,1.)
@test typeof(m) <: MaternCovarianceFunction
m = MaternCovarianceFunction(1,0.1,1.,2,1.)
@test typeof(m) <: MaternCovarianceFunction
m = MaternCovarianceFunction(1,0.1,1.,2.5,1.)
@test typeof(m) <: MaternCovarianceFunction
m = MaternCovarianceFunction(1,1,1,2.5,1.)
@test typeof(m) <: MaternCovarianceFunction
m = MaternCovarianceFunction(1,1,1.,2,1.)
@test typeof(m) <: MaternCovarianceFunction
m = MaternCovarianceFunction(1,1,1.,1,1)
@test typeof(m) <: MaternCovarianceFunction
m = MaternCovarianceFunction(2,0.1,1,2.5,1)
@test typeof(m) <: MaternCovarianceFunction
m = MaternCovarianceFunction(2,1//10,1,2,1)
@test typeof(m) <: MaternCovarianceFunction
m = MaternCovarianceFunction(2,1,1,2,1)
@test typeof(m) <: MaternCovarianceFunction

@test_throws ArgumentError MaternCovarianceFunction(1,-0.1,1.,2.5,1.)
@test_throws ArgumentError MaternCovarianceFunction(1,0.1,-1.,2.5,1.)
@test_throws ArgumentError MaternCovarianceFunction(1,0.1,1.,-2.5,1.)
@test_throws ArgumentError MaternCovarianceFunction(1,0.1,1.,2.5,0.5)
@test_throws ArgumentError MaternCovarianceFunction(0,0.1,1.,2.5,2)

verbose && println("done")

## SeparableMaternCovarianceFunction ##
verbose && print("testing separable matern covariance function...")

m = SeparableMaternCovarianceFunction(2,0.1,1.,2.5,1.)
@test typeof(m) <: SeparableMaternCovarianceFunction
m = SeparableMaternCovarianceFunction(2,1,1.,2.5,1.)
@test typeof(m) <: SeparableMaternCovarianceFunction
m = SeparableMaternCovarianceFunction(2,0.1,1,2.5,1.)
@test typeof(m) <: SeparableMaternCovarianceFunction
m = SeparableMaternCovarianceFunction(1,0.1,1.,2,1.)
@test typeof(m) <: SeparableMaternCovarianceFunction
m = SeparableMaternCovarianceFunction(1,0.1,1.,2.5,1.)
@test typeof(m) <: SeparableMaternCovarianceFunction
m = SeparableMaternCovarianceFunction(1,1,1,2.5,1.)
@test typeof(m) <: SeparableMaternCovarianceFunction
m = SeparableMaternCovarianceFunction(1,1,1.,2,1.)
@test typeof(m) <: SeparableMaternCovarianceFunction
m = SeparableMaternCovarianceFunction(1,1,1.,1,1)
@test typeof(m) <: SeparableMaternCovarianceFunction
m = SeparableMaternCovarianceFunction(2,0.1,1,2.5,1)
@test typeof(m) <: SeparableMaternCovarianceFunction
m = SeparableMaternCovarianceFunction(2,1//10,1,2,1)
@test typeof(m) <: SeparableMaternCovarianceFunction
m = SeparableMaternCovarianceFunction(2,1,1,2,1)
@test typeof(m) <: SeparableMaternCovarianceFunction

@test_throws ArgumentError SeparableMaternCovarianceFunction(1,-0.1,1.,2.5,1.)
@test_throws ArgumentError SeparableMaternCovarianceFunction(1,0.1,-1.,2.5,1.)
@test_throws ArgumentError SeparableMaternCovarianceFunction(1,0.1,1.,-2.5,1.)
@test_throws ArgumentError SeparableMaternCovarianceFunction(1,0.1,1.,2.5,0.5)
@test_throws ArgumentError SeparableMaternCovarianceFunction(0,0.1,1.,2.5,2)

verbose && println("done")

## ExponentialCovarianceFunction ##
verbose && print("testing exponential covariance function...")

m = ExponentialCovarianceFunction(2,0.1,1.,1.)
@test typeof(m) <: ExponentialCovarianceFunction
m = ExponentialCovarianceFunction(2,1,1.,1.)
@test typeof(m) <: ExponentialCovarianceFunction
m = ExponentialCovarianceFunction(2,0.1,1,1.)
@test typeof(m) <: ExponentialCovarianceFunction
m = ExponentialCovarianceFunction(1,0.1,1.,1.)
@test typeof(m) <: ExponentialCovarianceFunction
m = ExponentialCovarianceFunction(1,0.1,1.,1.)
@test typeof(m) <: ExponentialCovarianceFunction
m = ExponentialCovarianceFunction(1,1,1,1.)
@test typeof(m) <: ExponentialCovarianceFunction
m = ExponentialCovarianceFunction(1,1,1.,1.)
@test typeof(m) <: ExponentialCovarianceFunction
m = ExponentialCovarianceFunction(1,1,1.,1)
@test typeof(m) <: ExponentialCovarianceFunction
m = ExponentialCovarianceFunction(2,0.1,1,1)
@test typeof(m) <: ExponentialCovarianceFunction
m = ExponentialCovarianceFunction(2,1//10,1,1)
@test typeof(m) <: ExponentialCovarianceFunction
m = ExponentialCovarianceFunction(2,1,1,1)
@test typeof(m) <: ExponentialCovarianceFunction

@test_throws ArgumentError ExponentialCovarianceFunction(1,-0.1,1.,1.)
@test_throws ArgumentError ExponentialCovarianceFunction(1,0.1,-1,1.)
@test_throws ArgumentError ExponentialCovarianceFunction(1,0.1,1.,0.5)
@test_throws ArgumentError ExponentialCovarianceFunction(0,0.1,1.,2)

verbose && println("done")

## SeparableExponentialCovarianceFunction ##
verbose && print("testing separable exponential covariance function...")

m = SeparableExponentialCovarianceFunction(2,0.1,1.,1.)
@test typeof(m) <: SeparableExponentialCovarianceFunction
m = SeparableExponentialCovarianceFunction(2,1,1.,1.)
@test typeof(m) <: SeparableExponentialCovarianceFunction
m = SeparableExponentialCovarianceFunction(2,0.1,1,1.)
@test typeof(m) <: SeparableExponentialCovarianceFunction
m = SeparableExponentialCovarianceFunction(1,0.1,1.,1.)
@test typeof(m) <: SeparableExponentialCovarianceFunction
m = SeparableExponentialCovarianceFunction(1,0.1,1.,1.)
@test typeof(m) <: SeparableExponentialCovarianceFunction
m = SeparableExponentialCovarianceFunction(1,1,1,1.)
@test typeof(m) <: SeparableExponentialCovarianceFunction
m = SeparableExponentialCovarianceFunction(1,1,1.,1.)
@test typeof(m) <: SeparableExponentialCovarianceFunction
m = SeparableExponentialCovarianceFunction(1,1,1.,1)
@test typeof(m) <: SeparableExponentialCovarianceFunction
m = SeparableExponentialCovarianceFunction(2,0.1,1,1)
@test typeof(m) <: SeparableExponentialCovarianceFunction
m = SeparableExponentialCovarianceFunction(2,1//10,1,1)
@test typeof(m) <: SeparableExponentialCovarianceFunction
m = SeparableExponentialCovarianceFunction(2,1,1,1)
@test typeof(m) <: SeparableExponentialCovarianceFunction

@test_throws ArgumentError SeparableExponentialCovarianceFunction(1,-0.1,1.,1.)
@test_throws ArgumentError SeparableExponentialCovarianceFunction(1,0.1,-1.,1.)
@test_throws ArgumentError SeparableExponentialCovarianceFunction(1,0.1,1.,0.5)
@test_throws ArgumentError SeparableExponentialCovarianceFunction(0,0.1,1.,2)

verbose && println("done")
