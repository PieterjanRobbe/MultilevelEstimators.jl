# Number generator type
abstract NumberGenerator

abstract MCgenerator <: NumberGenerator

abstract QMCgenerator <: NumberGenerator

#
# Uniform random number generators
#

# uniform Monte Carlo sampler
type UniformMCgenerator{s} <: MCgenerator
  λ::Float64 # decay rate, 0.5 for Monte Carlo
  lb::Vector{Float64} # lower bounds for random variables
  ub::Vector{Float64} # upper bounds for random variables
end

# constructors
function UniformMCgenerator{N}(s::N)
  return UniformMCgenerator{s}(0.5,zeros(s),ones(s))
end

function UniformMCgenerator{N,T}(s::N,lb::Vector{T},ub::Vector{T})
  @assert length(ub) == s && length(lb) == s
  return UniformMCgenerator{s}(0.5,lb,ub)
end

# utilities
function getPoint{s,N}(generator::UniformMCgenerator{s},k::N)
  return generator.lb + (generator.ub - generator.lb).*rand(s)
end

nshifts{s}(generator::UniformMCgenerator{s}) = 1

reset{s}(generator::UniformMCgenerator{s}) = Void

# randomized QMC generator
type UniformQMCgenerator{s,q} <: QMCgenerator
  generator::RandWrapper
  λ::Float64 # decay rate, 1 for Rank-1 Lattice Rules
  lb::Vector{Float64} # lower bounds for random variables
  ub::Vector{Float64} # upper bounds for random variables
end

# constructors
function UniformQMCgenerator{N}(s::N, q::N)
  lat = LatSeq(s)
  randlat = RandWrapper(lat,q)
  return UniformQMCgenerator{s,q}(randlat,1.,zeros(s),ones(s))
end

function UniformQMCgenerator{N}(randlat::RandWrapper)
  return UniformQMCgenerator{s,q}(randlat,1.,zeros(s),ones(s))
end

function UniformQMCgenerator{N}(s::N, q::N,lb::Vector{T},ub::Vector{T})
  @assert length(ub) == s && length(lb) == s
  lat = LatSeq(s)
  randlat = RandWrapper(lat,q)
  return UniformQMCgenerator{s,q}(randlat,1.,lb,ub)
end

function UniformQMCgenerator{N}(randlat::RandWrapper,lb::Vector{T},ub::Vector{T})
  @assert length(ub) == s && length(lb) == s
  return UniformQMCgenerator{s,q}(randlat,1.,lb,ub)
end

# utilities
function getPoint{s,q,N}(generator::UniformQMCgenerator{s,q},k::N)
    return generator.lb + (generator.ub - generator.lb).*getPoint(generator.generator,k)
end

nshifts{s,q}(generator::UniformQMCgenerator{s,q}) = q

reset{s,q}(generator::UniformQMCgenerator{s,q})
  reset(generator.generator)
end

#
# Gaussian random number generators
#

# random number generator
type GaussianMCgenerator{s} <: MCgenerator
  λ::Float64
end

# constructors
function GaussianMCgenerator{N}(s::N)
  return GaussianMCgenerator{s,1}(0.5)
end

# utilities
getPoint{s,N}(generator::GaussianMCgenerator{s},k::N) = randn(s)

nshifts{s}(generator::GaussianMCgenerator{s}) = 1

reset{s}(generator::GaussianMCgenerator{s}) = Void

# randomised lattice rule generator
type GaussianQMCgenerator{s,q} <: QMCgenerator
  generator::RandWrapper
  λ::Float64
end

# constructors
function GaussianQMCgenerator{N}(s::N, q::N)
  lat = LatSeq(s)
  randlat = RandWrapper(lat,q)
  return GaussianQMCgenerator{s,q}(randlat,1.)
end

function GaussianQMCgenerator{N}(randlat::RandWrapper)
  return GaussianQMCgenerator{s,q}(randlat,1.)
end

# utilities
function getPoint{s,q,N}(generator::GaussianQMCgenerator{s,q},k::N)
    return getPoint(generator.generator,k)
end

nshifts{s,q}(generator::GaussianQMCgenerator{s,q}) = q

reset{s,q}(generator::GaussianQMCgenerator{s,q})
  reset(generator.generator)
end