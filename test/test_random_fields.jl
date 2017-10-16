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
m = MaternCovarianceFunction(2,0.1,1.,2.5,1.,is_separable=false)
@test typeof(m) <: MaternCovarianceFunction
m = ExponentialCovarianceFunction(2,0.1,1.,1.,is_separable=false)
@test typeof(m) <: MaternCovarianceFunction

@test_throws ArgumentError MaternCovarianceFunction(1,-0.1,1.,2.5,1.)
@test_throws ArgumentError MaternCovarianceFunction(1,0.1,-1.,2.5,1.)
@test_throws ArgumentError MaternCovarianceFunction(1,0.1,1.,-2.5,1.)
@test_throws ArgumentError MaternCovarianceFunction(1,0.1,1.,2.5,0.5)
@test_throws ArgumentError MaternCovarianceFunction(0,0.1,1.,2.5,2)

verbose && println("done")
