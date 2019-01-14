## distribution.jl : unit testing for distribution.jl
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
@sprintf "%s" distr1
for xᵢ in 0.1:0.1:0.9
    @test transform(Uniform(), xᵢ) == xᵢ
end

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
@sprintf "%s" distr1
y = [-1.2815515655446004, -0.8416212335729143, -0.5244005127080409, -0.2533471031357997, 0.0,
      0.2533471031357997,  0.5244005127080407,  0.8416212335729143,  1.2815515655446004]
for (xᵢ, yᵢ) in Base.Iterators.zip(0.1:0.1:0.9, y)
    @test transform(Normal(), xᵢ) ≈ yᵢ
end

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
@sprintf "%s" distr1
y = [ -1.184032466693905, -0.7938201191289552, -0.4984028824219386, -0.2415871851410768, 0.0,
      0.2415871851410768,  0.4984028824219383,  0.7938201191289559,  1.1840324666939055] 
for (xᵢ, yᵢ) in Base.Iterators.zip(0.1:0.1:0.9, y)
    @test transform(TruncatedNormal(), xᵢ) ≈ yᵢ
end

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
@sprintf "%s" distr1
y = [  0.4590436050264207, 0.6680472308365775,  0.8446004309005917,  1.0107676525947897, 1.1774100225154747,
       1.3537287260556712, 1.5517556536555206,  1.7941225779941017,  2.1459660262893476]
for (xᵢ, yᵢ) in Base.Iterators.zip(0.1:0.1:0.9, y)
    @test transform(Weibull(), xᵢ) ≈ yᵢ
end

end
