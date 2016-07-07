# Number generator type
abstract NumberGenerator

abstract MCgenerator <: NumberGenerator

abstract QMCgenerator <: NumberGenerator

# Uniform random number generators

# random number generator
type UniformMCgenerator{s,q,d} <: MCgenerator
  λ::Float64
end

function UniformMCgenerator{N}(d::N, s::N)
  return UniformMCgenerator{s,1,d}(0.5)
end

next{s,q,d}(generator::UniformMCgenerator{s,q,d},index::Index{d}) = rand(s)

# randomised lattice rule generator
type UniformQMCgenerator{s,q,d} <: QMCgenerator
  generators::Dict{Index{d},RandLatSeq{s,q}}
  λ::Float64
end

function UniformQMCgenerator{N}(d::N, s::N, q::N)
  return UniformQMCgenerator{s,q,d}(Dict{Index{d},RandLatSeq{s,q}}(),1.)
end

function next{s,q,d}(qmc::UniformQMCgenerator{s,q,d},index::Index{d})
  if haskey(qmc.generators,index)
    next!(qmc.generators[index])
  else
    qmc.generators[index] = RandLatSeq(s,q)
    next!(qmc.generators[index])
  end
end

nshifts{s,q,d}(generator::UniformMCgenerator{s,q,d}) = q
nshifts{s,q,d}(generator::UniformQMCgenerator{s,q,d}) = q

function reset!{s,q,d}(generator::UniformMCgenerator{s,q,d})
end

function reset!{s,q,d}(generator::UniformQMCgenerator{s,q,d})
  for gen in values(generator.generators)
    QMC.reset!(gen)
  end
end

# Gaussian random number generators

# random number generator
type GaussianMCgenerator{s,q,d} <: MCgenerator
  λ::Float64
end

function GaussianMCgenerator{N}(d::N, s::N)
  return GaussianMCgenerator{s,1,d}(0.5)
end

next{s,q,d}(generator::GaussianMCgenerator{s,q,d},index::Index{d}) = randn(s)

# randomised lattice rule generator
type GaussianQMCgenerator{s,q,d} <: QMCgenerator
  generators::Dict{Index{d},RandLatSeq{s,q}}
  λ::Float64
end

function GaussianQMCgenerator{N}(d::N, s::N, q::N)
  return GaussianQMCgenerator{s,q,d}(Dict{Index{d},RandLatSeq{s,q}}(),1.)
end

function next{s,q,d}(qmc::GaussianQMCgenerator{s,q,d},index::Index{d})
if haskey(qmc.generators,index)
    sqrt(2)*erfinv(2*next!(qmc.generators[index])-1)
  else
    qmc.generators[index] = RandLatSeq(s,q)
    sqrt(2)*erfinv(2*next!(qmc.generators[index])-1)
  end
end

nshifts{s,q,d}(generator::GaussianMCgenerator{s,q,d}) = q
nshifts{s,q,d}(generator::GaussianQMCgenerator{s,q,d}) = q

function reset!{s,q,d}(generator::GaussianMCgenerator{s,q,d})
end

function reset!{s,q,d}(generator::GaussianQMCgenerator{s,q,d})
  for gen in values(generator.generators)
    QMC.reset!(gen)
  end
end