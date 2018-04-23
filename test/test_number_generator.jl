# test_number_generator.jl : test NumberGenerator

@testset "NumberGenerator              " begin

# UniformMCGenerator
u = UniformMCGenerator(100)
@test typeof(u) <: MultilevelEstimators.NumberGenerator
u = UniformMCGenerator(-1*ones(100),ones(100))
@test typeof(u) <: MultilevelEstimators.NumberGenerator

# UniformQMCGenerator
u = UniformQMCGenerator(100,16)
@test typeof(u) <: MultilevelEstimators.NumberGenerator
u = UniformQMCGenerator(16,-1*ones(100),ones(100))
@test typeof(u) <: MultilevelEstimators.NumberGenerator
lat = LatSeq(12)
u = UniformQMCGenerator(lat,16,-1*ones(12),ones(12))
@test typeof(u) <: MultilevelEstimators.NumberGenerator

# NormalMCGenerator
n = NormalMCGenerator(250)
@test typeof(n) <: MultilevelEstimators.NumberGenerator
n = NormalMCGenerator(-3*ones(101),2*ones(101))
@test typeof(n) <: MultilevelEstimators.NumberGenerator

# NormalQMCGenerator
n = NormalQMCGenerator(100,16)
@test typeof(n) <: MultilevelEstimators.NumberGenerator
n = NormalQMCGenerator(16,zeros(100),ones(100))
@test typeof(n) <: MultilevelEstimators.NumberGenerator
lat = LatSeq(100)
n = NormalQMCGenerator(lat,64,-1*ones(100),ones(100))
@test typeof(n) <: MultilevelEstimators.NumberGenerator

# TruncatedNormalMCGenerator
t = TruncatedNormalMCGenerator(25)
@test typeof(t) <: MultilevelEstimators.NumberGenerator
t = TruncatedNormalMCGenerator(-3*ones(101),2*ones(101),-5*ones(101),-1*ones(101))
@test typeof(t) <: MultilevelEstimators.NumberGenerator

# TruncatedNormalQMCGenerator
t = TruncatedNormalQMCGenerator(10,16)
@test typeof(t) <: MultilevelEstimators.NumberGenerator
t = TruncatedNormalQMCGenerator(16,zeros(100),ones(100),-2*ones(100),2*ones(100))
@test typeof(t) <: MultilevelEstimators.NumberGenerator
lat = LatSeq(3500)
t = TruncatedNormalQMCGenerator(lat,16,-1*ones(3500),ones(3500),-4*ones(3500),ones(3500))
@test typeof(t) <: MultilevelEstimators.NumberGenerator

# argument error handling
@test_throws BoundsError UniformMCGenerator(0) 
@test_throws BoundsError UniformQMCGenerator(10,0) 
@test_throws ArgumentError UniformMCGenerator(zeros(300)) 
@test_throws ArgumentError UniformMCGenerator(zeros(300),ones(300),zeros(300)) 
@test_throws ArgumentError NormalMCGenerator(zeros(300),ones(300),zeros(300)) 
@test_throws ArgumentError NormalMCGenerator(-ones(300),-ones(300)) 
@test_throws ArgumentError TruncatedNormalMCGenerator(zeros(300),ones(300),zeros(300)) 
@test_throws ArgumentError TruncatedNormalMCGenerator(zeros(300),ones(300)) 
@test_throws ArgumentError TruncatedNormalMCGenerator(zeros(300),-ones(300),zeros(300),ones(300)) 
@test_throws ArgumentError UniformMCGenerator(ones(300),zeros(300)) 
@test_throws DimensionMismatch UniformMCGenerator(zeros(300),ones(301)) 

end
