## distributions.jl : unit testing for distributions.jl
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

@testset "Uniform                      " begin

distr1 = Uniform()
distr2 = Uniform(0.,1.)
distr3 = Uniform(0,1.)
@test typeof(distr1) <: Uniform
@test typeof(distr2) <: Uniform
@test typeof(distr3) <: Uniform
@test_throws ArgumentError Uniform(Inf, 1.)
@test_throws ArgumentError Uniform(0., NaN)
@test_throws ArgumentError Uniform(1., 0.)
str = @sprintf "%s" distr1

end

@testset "Normal                       " begin

distr1 = Normal()
distr2 = Normal(0.,1.)
distr3 = Normal(0,1.)
@test typeof(distr1) <: Normal
@test typeof(distr2) <: Normal
@test typeof(distr3) <: Normal
@test_throws ArgumentError Normal(Inf, 1.)
@test_throws ArgumentError Normal(0., NaN)
@test_throws ArgumentError Normal(0., -1.)
str = @sprintf "%s" distr1

end

@testset "TruncatedNormal              " begin

distr1 = TruncatedNormal()
distr2 = TruncatedNormal(0.,1.)
distr3 = TruncatedNormal(0.,1.,-1.,1.)
distr4 = TruncatedNormal(0.,1.,-2//3,10)
@test typeof(distr1) <: TruncatedNormal
@test typeof(distr2) <: TruncatedNormal
@test typeof(distr3) <: TruncatedNormal
@test typeof(distr4) <: TruncatedNormal
@test_throws ArgumentError TruncatedNormal(Inf, 1.)
@test_throws ArgumentError TruncatedNormal(0., NaN)
@test_throws ArgumentError TruncatedNormal(0., 1., Inf, 2.)
@test_throws ArgumentError TruncatedNormal(0., 1., -2., Inf)
@test_throws ArgumentError TruncatedNormal(0., -1.)
@test_throws ArgumentError TruncatedNormal(1., 0.)
str = @sprintf "%s" distr1

end

@testset "Weibull                      " begin

distr1 = Weibull()
distr2 = Weibull(2.,1.)
distr3 = Weibull(2,1.)
@test typeof(distr1) <: Weibull
@test typeof(distr2) <: Weibull
@test typeof(distr3) <: Weibull
@test_throws ArgumentError Weibull(Inf, 1.)
@test_throws ArgumentError Weibull(0., NaN)
@test_throws ArgumentError Weibull(1., -1)
@test_throws ArgumentError Weibull(-1., 1)
str = @sprintf "%s" distr1

end
