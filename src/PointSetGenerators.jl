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
function getPoint{s,N<:Integer}(generator::UniformMCgenerator{s},k::N)
  return generator.lb + (generator.ub - generator.lb).*rand(s)
end

nshifts{s}(generator::UniformMCgenerator{s}) = 1

reset{s}(generator::UniformMCgenerator{s}) = Void

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
  return UniformQMCgenerator{s,q,typeof(randlat),Float64,Vector{Float64}}(randlat,1.,zeros(s),ones(s))
end

function UniformQMCgenerator{N<:Integer,T<:AbstractFloat}(s::N, q::N,lb::Vector{T},ub::Vector{T})
  @assert length(ub) == s && length(lb) == s
  lat = LatSeq(s)
  randlat = RandWrapper(lat,q)
  return UniformQMCgenerator{s,q,typeof(randlat),T,Vector{T}}(randlat,1.,lb,ub)
end

function UniformQMCgenerator{T<:AbstractFloat}(randlat::RandWrapper,lb::Vector{T},ub::Vector{T})
  @assert length(ub) == s && length(lb) == s
  return UniformQMCgenerator{s,q,typeof(randlat),T,Vector{T}}(randlat,1.,lb,ub)
end

# utilities
function getPoint{s,q,N<:Integer}(generator::UniformQMCgenerator{s,q},k::N)
    return generator.lb + (generator.ub - generator.lb).*getPoint(generator.generator,k)
end

nshifts{s,q}(generator::UniformQMCgenerator{s,q}) = q

reset(generator::UniformQMCgenerator) = reset(generator.generator)

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
getPoint{s,N<:Integer}(generator::GaussianMCgenerator{s},k::N) = randn(s)

nshifts{s}(generator::GaussianMCgenerator{s}) = 1

reset{s}(generator::GaussianMCgenerator{s}) = Void

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
function getPoint{s,q,N}(generator::GaussianQMCgenerator{s,q},k::N)
    return sqrt(2)*erfinv(2*getPoint(generator.generator,k)-1)
end

nshifts{s,q}(generator::GaussianQMCgenerator{s,q}) = q

reset(generator::GaussianQMCgenerator) = reset(generator.generator)

ndims{s}(generator::NumberGenerator{s}) = s