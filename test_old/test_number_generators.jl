## test_number_generators.jl : tests for number_generators.jl

## uniform MC generator ##
verbose && print("testing uniform MC generator...")

umc = UniformMCgenerator(200)
@test typeof(umc) <: UniformMCgenerator
umc = UniformMCgenerator(20,-ones(20),ones(20))
@test typeof(umc) <: UniformMCgenerator
@test_throws ArgumentError UniformMCgenerator(-10)
@test_throws ArgumentError UniformMCgenerator(10,-ones(9),ones(10))
@test_throws ArgumentError UniformMCgenerator(10,-ones(10),ones(9))
@test_throws ArgumentError UniformMCgenerator(10,ones(10),-ones(10))

verbose && println("done")

## uniform QMC generator ##
verbose && print("testing uniform QMC generator...")

uqmc = UniformQMCgenerator(200,16)
@test typeof(uqmc) <: UniformQMCgenerator
uqmc = UniformQMCgenerator(20,8,-ones(20),ones(20))
@test typeof(uqmc) <: UniformQMCgenerator
lat = LatSeq(16)
randlat = RandWrapper(lat,8)
uqmc = UniformQMCgenerator(randlat)
@test typeof(uqmc) <: UniformQMCgenerator
dig = DigSeq(16)
randdig = RandWrapper(dig,8)
uqmc = UniformQMCgenerator(randdig)
@test typeof(uqmc) <: UniformQMCgenerator
s = 16; q = 8
lat = LatSeq(s)
randlat = RandWrapper(lat,q)
uqmc = UniformQMCgenerator(randlat,-ones(s),ones(s))
@test typeof(uqmc) <: UniformQMCgenerator
@test_throws ArgumentError UniformQMCgenerator(-10,4)
@test_throws ArgumentError UniformQMCgenerator(10,-4)
@test_throws ArgumentError UniformQMCgenerator(10,4,-ones(9),ones(10))
@test_throws ArgumentError UniformQMCgenerator(10,4,-ones(10),ones(9))
@test_throws ArgumentError UniformQMCgenerator(10,4,ones(10),-ones(10))
@test_throws ArgumentError UniformQMCgenerator(randlat,-ones(s-1),ones(s))
@test_throws ArgumentError UniformQMCgenerator(randlat,-ones(s),ones(s-1))
@test_throws ArgumentError UniformQMCgenerator(randlat,ones(s),-ones(s))

verbose && println("done")

## Gaussian MC generator ##
verbose && print("testing Gaussian MC generator...")

gmc = GaussianMCgenerator(200)
@test typeof(gmc) <: GaussianMCgenerator
@test_throws ArgumentError GaussianMCgenerator(-10)

verbose && println("done")

## Gaussian QMC generator ##
verbose && print("testing Gaussian QMC generator...")

gqmc = GaussianQMCgenerator(800,16)
@test typeof(gqmc) <: GaussianQMCgenerator
lat = LatSeq(16)
randlat = RandWrapper(lat,8)
gqmc = GaussianQMCgenerator(randlat)
@test typeof(gqmc) <: GaussianQMCgenerator
@test_throws ArgumentError GaussianQMCgenerator(-800,16)
@test_throws ArgumentError GaussianQMCgenerator(800,-16)

verbose && println("done")

