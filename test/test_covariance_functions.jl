# TODO \'e
## test_covariance_functions.jl : tests for covariance_functions.jl

## MaternCovarianceFunction ##
verbose && print("testing Mat\'ern covariance function...")

m = MaternCovarianceFunction(2,0.1,2.5,1.)
@test typeof(m) <: CovarianceFunction
m = MaternCovarianceFunction(2,1,2.5,1.)
@test typeof(m) <: CovarianceFunction
m = MaternCovarianceFunction(2,0.1,2.5,1.)
@test typeof(m) <: CovarianceFunction
m = MaternCovarianceFunction(1,0.1,2,1.)
@test typeof(m) <: CovarianceFunction
m = MaternCovarianceFunction(1,0.1,2.5,1.)
@test typeof(m) <: CovarianceFunction
m = MaternCovarianceFunction(1,1,2.5,1.)
@test typeof(m) <: CovarianceFunction
m = MaternCovarianceFunction(1,1,2,1.)
@test typeof(m) <: CovarianceFunction
m = MaternCovarianceFunction(1,1,1,1)
@test typeof(m) <: CovarianceFunction
m = MaternCovarianceFunction(2,0.1,2.5,1)
@test typeof(m) <: CovarianceFunction
m = MaternCovarianceFunction(2,1//10,2,1)
@test typeof(m) <: CovarianceFunction
m = MaternCovarianceFunction(2,1,2,1)
@test typeof(m) <: CovarianceFunction

@test_throws ArgumentError MaternCovarianceFunction(1,-0.1,2.5,1.)
@test_throws ArgumentError MaternCovarianceFunction(1,0.1,-2.5,1.)
@test_throws ArgumentError MaternCovarianceFunction(1,0.1,2.5,0.5)
@test_throws ArgumentError MaternCovarianceFunction(0,0.1,2.5,2)
@test_throws MethodError MaternCovarianceFunction(0.1,0.1,2.5,2)

verbose && println("done")

## SeparableMaternCovarianceFunction ##
verbose && print("testing separable Mat\'ern covariance function...")

m = SeparableMaternCovarianceFunction(2,0.1,2.5,1.)
@test typeof(m) <: CovarianceFunction
m = SeparableMaternCovarianceFunction(2,1,2.5,1.)
@test typeof(m) <: CovarianceFunction
m = SeparableMaternCovarianceFunction(2,0.1,2.5,1.)
@test typeof(m) <: CovarianceFunction
m = SeparableMaternCovarianceFunction(1,0.1,2,1.)
@test typeof(m) <: CovarianceFunction
m = SeparableMaternCovarianceFunction(1,0.1,2.5,1.)
@test typeof(m) <: CovarianceFunction
m = SeparableMaternCovarianceFunction(1,1,2.5,1.)
@test typeof(m) <: CovarianceFunction
m = SeparableMaternCovarianceFunction(1,1,2,1.)
@test typeof(m) <: CovarianceFunction
m = SeparableMaternCovarianceFunction(1,1,1,1)
@test typeof(m) <: CovarianceFunction
m = SeparableMaternCovarianceFunction(2,0.1,2.5,1)
@test typeof(m) <: CovarianceFunction
m = SeparableMaternCovarianceFunction(2,1//10,2,1)
@test typeof(m) <: CovarianceFunction
m = SeparableMaternCovarianceFunction(2,1,2,1)
@test typeof(m) <: CovarianceFunction

@test_throws ArgumentError SeparableMaternCovarianceFunction(1,-0.1,2.5,1.)
@test_throws ArgumentError SeparableMaternCovarianceFunction(1,0.1,-2.5,1.)
@test_throws ArgumentError SeparableMaternCovarianceFunction(1,0.1,2.5,0.5)
@test_throws ArgumentError SeparableMaternCovarianceFunction(0,0.1,1.,2)
@test_throws MethodError SeparableMaternCovarianceFunction(0.1,0.1,2.5,2)

verbose && println("done")

## ExponentialCovarianceFunction ##
verbose && print("testing exponential covariance function...")

e = ExponentialCovarianceFunction(2,0.1,1.)
@test typeof(e) <: CovarianceFunction
e = ExponentialCovarianceFunction(2,1,1.)
@test typeof(e) <: CovarianceFunction
e = ExponentialCovarianceFunction(2,0.1,1.)
@test typeof(e) <: CovarianceFunction
e = ExponentialCovarianceFunction(1,0.1,1.)
@test typeof(e) <: CovarianceFunction
e = ExponentialCovarianceFunction(1,0.1,1.)
@test typeof(e) <: CovarianceFunction
e = ExponentialCovarianceFunction(1,1,1.)
@test typeof(e) <: CovarianceFunction
e = ExponentialCovarianceFunction(1,1,1.)
@test typeof(e) <: CovarianceFunction
e = ExponentialCovarianceFunction(1,1,1)
@test typeof(e) <: CovarianceFunction
e = ExponentialCovarianceFunction(2,0.1,1)
@test typeof(e) <: CovarianceFunction
e = ExponentialCovarianceFunction(2,1//10,1)
@test typeof(e) <: CovarianceFunction
e = ExponentialCovarianceFunction(2,1,1)
@test typeof(e) <: CovarianceFunction

@test_throws ArgumentError ExponentialCovarianceFunction(1,-0.1,1.)
@test_throws ArgumentError ExponentialCovarianceFunction(1,0.1,0.5)
@test_throws ArgumentError ExponentialCovarianceFunction(0,0.1,2)
@test_throws MethodError ExponentialCovarianceFunction(0.1,0.1,2)

verbose && println("done")

## SeparableExponentialCovarianceFunction ##
verbose && print("testing separable exponential covariance function...")

e = SeparableExponentialCovarianceFunction(2,0.1,1.)
@test typeof(e) <: CovarianceFunction
e = SeparableExponentialCovarianceFunction(2,1,1.)
@test typeof(e) <: CovarianceFunction
e = SeparableExponentialCovarianceFunction(2,0.1,1.)
@test typeof(e) <: CovarianceFunction
e = SeparableExponentialCovarianceFunction(1,0.1,1.)
@test typeof(e) <: CovarianceFunction
e = SeparableExponentialCovarianceFunction(1,0.1,1.)
@test typeof(e) <: CovarianceFunction
e = SeparableExponentialCovarianceFunction(1,1,1)
@test typeof(e) <: CovarianceFunction
e = SeparableExponentialCovarianceFunction(1,1,1.)
@test typeof(e) <: CovarianceFunction
e = SeparableExponentialCovarianceFunction(1,1,1)
@test typeof(e) <: CovarianceFunction
e = SeparableExponentialCovarianceFunction(2,0.1,1)
@test typeof(e) <: CovarianceFunction
e = SeparableExponentialCovarianceFunction(2,1//10,1)
@test typeof(e) <: CovarianceFunction
e = SeparableExponentialCovarianceFunction(2,1,1)
@test typeof(e) <: CovarianceFunction

@test_throws ArgumentError SeparableExponentialCovarianceFunction(1,-0.1,1.)
@test_throws ArgumentError SeparableExponentialCovarianceFunction(1,0.1,0.5)
@test_throws ArgumentError SeparableExponentialCovarianceFunction(0,0.1,2)
@test_throws MethodError SeparableExponentialCovarianceFunction(0.1,0.1,2)

verbose && println("done")

## SquaredExponentialCovarianceFunction ##
verbose && print("testing squared exponential covariance function...")

e = SquaredExponentialCovarianceFunction(2,0.1,1.)
@test typeof(e) <: CovarianceFunction
e = SquaredExponentialCovarianceFunction(2,1,1.)
@test typeof(e) <: CovarianceFunction
e = SquaredExponentialCovarianceFunction(2,0.1,1.)
@test typeof(e) <: CovarianceFunction
e = SquaredExponentialCovarianceFunction(1,0.1,1.)
@test typeof(e) <: CovarianceFunction
e = SquaredExponentialCovarianceFunction(1,0.1,1.)
@test typeof(e) <: CovarianceFunction
e = SquaredExponentialCovarianceFunction(1,1,1.)
@test typeof(e) <: CovarianceFunction
e = SquaredExponentialCovarianceFunction(1,1,1.)
@test typeof(e) <: CovarianceFunction
e = SquaredExponentialCovarianceFunction(1,1,1)
@test typeof(e) <: CovarianceFunction
e = SquaredExponentialCovarianceFunction(2,0.1,1)
@test typeof(e) <: CovarianceFunction
e = SquaredExponentialCovarianceFunction(2,1//10,1)
@test typeof(e) <: CovarianceFunction
e = SquaredExponentialCovarianceFunction(2,1,1)
@test typeof(e) <: CovarianceFunction

@test_throws ArgumentError SquaredExponentialCovarianceFunction(1,-0.1,1.)
@test_throws ArgumentError SquaredExponentialCovarianceFunction(1,0.1,0.5)
@test_throws ArgumentError SquaredExponentialCovarianceFunction(0,0.1,2)
@test_throws MethodError SquaredExponentialCovarianceFunction(0.1,0.1,2)

verbose && println("done")

## SeparableSquaredExponentialCovarianceFunction ##
verbose && print("testing separable squared exponential covariance function...")

e = SeparableSquaredExponentialCovarianceFunction(2,0.1,1.)
@test typeof(e) <: CovarianceFunction
e = SeparableSquaredExponentialCovarianceFunction(2,1,1.)
@test typeof(e) <: CovarianceFunction
e = SeparableSquaredExponentialCovarianceFunction(2,0.1,1.)
@test typeof(e) <: CovarianceFunction
e = SeparableSquaredExponentialCovarianceFunction(1,0.1,1.)
@test typeof(e) <: CovarianceFunction
e = SeparableSquaredExponentialCovarianceFunction(1,0.1,1.)
@test typeof(e) <: CovarianceFunction
e = SeparableSquaredExponentialCovarianceFunction(1,1,1)
@test typeof(e) <: CovarianceFunction
e = SeparableSquaredExponentialCovarianceFunction(1,1,1.)
@test typeof(e) <: CovarianceFunction
e = SeparableSquaredExponentialCovarianceFunction(1,1,1)
@test typeof(e) <: CovarianceFunction
e = SeparableSquaredExponentialCovarianceFunction(2,0.1,1)
@test typeof(e) <: CovarianceFunction
e = SeparableSquaredExponentialCovarianceFunction(2,1//10,1)
@test typeof(e) <: CovarianceFunction
e = SeparableSquaredExponentialCovarianceFunction(2,1,1)
@test typeof(e) <: CovarianceFunction

@test_throws ArgumentError SeparableSquaredExponentialCovarianceFunction(1,-0.1,1.)
@test_throws ArgumentError SeparableSquaredExponentialCovarianceFunction(1,0.1,0.5)
@test_throws ArgumentError SeparableSquaredExponentialCovarianceFunction(0,0.1,2)
@test_throws MethodError SeparableSquaredExponentialCovarianceFunction(0.1,0.1,2)

verbose && println("done")
