# Number generator type
abstract NumberGenerator{s}

abstract MCgenerator{s} <: NumberGenerator{s}

abstract QMCgenerator{s} <: NumberGenerator{s}

#
# Uniform random number generators
#

# uniform Monte Carlo sampler
type UniformMCgenerator{s,T<:AbstractFloat,V<:AbstractVector} <: MCgenerator{s}
  λ::T # decay rate, 0.5 for Monte Carlo
  lb::V # lower bounds for random variables
  ub::V # upper bounds for random variables
end

# constructors
function UniformMCgenerator{N<:Integer}(s::N)
  return UniformMCgenerator{s,Float64,Vector{Float64}}(0.5,zeros(s),ones(s))
end

function UniformMCgenerator{N<:Integer,T<:AbstractFloat}(s::N,lb::Vector{T},ub::Vector{T})
  @assert length(ub) == s && length(lb) == s
  return UniformMCgenerator{s,T,Vector{T}}(convert(T,0.5),lb,ub)
end

# utilities
nshifts{s}(generator::UniformMCgenerator{s}) = 1
reset{s}(generator::UniformMCgenerator{s}) = Void

# methods
isValid{s}(U::UniformMCgenerator{s}) = G.λ == 0.5 && length(lb) == s && length(ub) == s

function show{s}(io::IO,U::UniformMCgenerator{s})
  @assert isValid(U)
  str = "$(ndims(G))-dimensional uniform Monte Carlo sampler with λ= $(U.λ) and \n"
  str *= " - lower bound = $(U.lb)\n"
  str *= " - upper bound = $(U.lb)\n"
  print(io,str)
end

function getPoint{s,N<:Integer}(generator::UniformMCgenerator{s},k::N)
  return generator.lb .+ (generator.ub - generator.lb).*rand(s)
end

# randomized QMC generator
type UniformQMCgenerator{s,q,R<:RandWrapper,T<:AbstractFloat,V<:AbstractVector} <: QMCgenerator{s}
  generator::R
  λ::T # decay rate, 1 for rank-1 lattice rules
  lb::V # lower bounds for random variables
  ub::V # upper bounds for random variables
end

# constructors
function UniformQMCgenerator{N<:Integer}(s::N, q::N)
  lat = LatSeq(s)
  randlat = RandWrapper(lat,q)
  return UniformQMCgenerator{s,q,typeof(randlat),Float64,Vector{Float64}}(randlat,1.,zeros(s),ones(s))
end

function UniformQMCgenerator(randlat::RandWrapper)
  return UniformQMCgenerator{ndims(s),nshifts(q),typeof(randlat),Float64,Vector{Float64}}(randlat,1.,zeros(s),ones(s))
end

function UniformQMCgenerator{N<:Integer,T<:AbstractFloat}(s::N, q::N,lb::Vector{T},ub::Vector{T})
  @assert length(ub) == s && length(lb) == s
  lat = LatSeq(s)
  randlat = RandWrapper(lat,q)
  return UniformQMCgenerator{s,q,typeof(randlat),T,Vector{T}}(randlat,1.,lb,ub)
end

function UniformQMCgenerator{T<:AbstractFloat}(randlat::RandWrapper,lb::Vector{T},ub::Vector{T})
  @assert length(ub) == ndims(randlat) && length(lb) == ndims(randlat)
  return UniformQMCgenerator{ndims(randlat),q = nshifts(randlat),typeof(randlat),T,Vector{T}}(randlat,1.,lb,ub)
end

function UniformQMCgenerator{T<:AbstractFloat}(randlat::RandWrapper,λ::T,lb::Vector{T},ub::Vector{T})
  @assert length(ub) == s && length(lb) == s
  return UniformQMCgenerator{ndims(randlat),q = nshifts(randlat),typeof(randlat),T,Vector{T}}(randlat,λ,lb,ub)
end

# utilities
nshifts{s,q}(generator::UniformQMCgenerator{s,q}) = q
reset(generator::UniformQMCgenerator) = reset(generator.generator)

# methods
isValid{s}(U::UniformQMCgenerator{s}) = G.λ >= 0.5 && length(lb) == s && length(ub) == s

function show{s,q}(io::IO,U::UniformQMCgenerator{s,q})
  @assert isValid(U)
  str = "$(nshifts(U)) x $(ndims(U))-dimensional uniform Quasi-Monte Carlo sampler with λ= $(U.λ) and \n"
  str *= " - lower bound = $(U.lb)\n"
  str *= " - upper bound = $(U.lb)\n"
  print(io,str)
end

function getPoint{s,q,N<:Integer}(generator::UniformQMCgenerator{s,q},k::N)
    return generator.lb .+ (generator.ub - generator.lb).*getPoint(generator.generator,k)
end

#
# Gaussian random number generators
#

# random number generator
type GaussianMCgenerator{s,T<:AbstractFloat} <: MCgenerator{s}
  λ::T
end

# constructors
function GaussianMCgenerator{N<:Integer}(s::N)
  return GaussianMCgenerator{s,Float64}(0.5)
end

# utilities
nshifts{s}(generator::GaussianMCgenerator{s}) = 1
reset{s}(generator::GaussianMCgenerator{s}) = Void

# methods
isValid(G::GaussianMCgenerator) = G.λ == 0.5

function show{s}(io::IO,G::GaussianMCgenerator{s})
  @assert isValid(G)
  str = "$(ndims(G))-dimensional Gaussian Monte Carlo sampler with λ= $(G.λ)"
  print(io,str)
end

getPoint{s,N<:Integer}(generator::GaussianMCgenerator{s},k::N) = randn(s)

# randomised lattice rule generator
type GaussianQMCgenerator{s,q,R<:RandWrapper,T<:AbstractFloat} <: QMCgenerator{s}
  generator::R
  λ::T
end

# constructors
function GaussianQMCgenerator{N<:Integer}(s::N, q::N)
  lat = LatSeq(s)
  randlat = RandWrapper(lat,q)
  return GaussianQMCgenerator{s,q,typeof(randlat),Float64}(randlat,1.)
end

function GaussianQMCgenerator(randlat::RandWrapper)
  return GaussianQMCgenerator{s,q,typeof(randlat),Float64}(randlat,1.)
end

# utilities
nshifts{s,q}(generator::GaussianQMCgenerator{s,q}) = q
ndims{s}(generator::NumberGenerator{s}) = s

# methods
isValid(G::GaussianQMCgenerator) = G.λ >= 0.5

function show{s,q}(io::IO,G::GaussianQMCgenerator{s,q})
  @assert isValid(G)
  str = "$(nshifts(G)) x $(ndims(G))-dimensional Gaussian Quasi-Monte Carlo sampler with λ= $(G.λ)"
  print(io,str)
end

function getPoint{s,q,N}(generator::GaussianQMCgenerator{s,q},k::N)
    return sqrt(2)*erfinv(2*getPoint(generator.generator,k)-1)
end

# reset(generator::GaussianQMCgenerator) = reset(generator.generator)














