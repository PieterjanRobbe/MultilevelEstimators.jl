# settings type
type Settings{d,V<:AbstractVector,N<:Integer,T<:AbstractFloat,I<:IndexSet,G<:NumberGenerator,F<:GaussianFieldSampler,D<:Dict}
  m0::V # coarsest mesh size
  maxK::N # max index number
  gamma::T # solver complexity
  splitting::T # desired MSE splitting
  Z::N # number of QOI's to determine
  indexset::I # type of index set
  numberGenerator::G # type of point generator
  gaussianFieldSampler::F # type of Gaussian field sampler
  sampleFunction::Function # quantity of interest
  showInfo::Bool # boolean to decide if something must be printed to screen
  useTime::Bool # boolean to indicate if true simulation time must be used as cost measure
  generatorState::D # save states of number generator
end

# utilities
inttype{d,V,N}(settings::Settings{d,V,N}) = N

ndims{d}(settings::Settings{d}) = d 

# outer constructor
function Settings{I,G,V,N,T,F}(
  indexset::I,
  numberGenerator::G,
  sampleFunction::Function;
  m0::V = 4*ones(Int64,1),
  maxK::N = 6,
  gamma::T = 1.4,
  splitting::T = 0.5,
  Z::N = 1, 
  gaussianFieldSampler::F = EmptySampler(),
  showInfo::Bool = true,
  useTime::Bool = false )
  d = ndims(indexset)
  while length(m0) < d # ensure that m0 has the correct dimension, since we don't know d on method call
    push!(m0,4)
  end
  generatorState = Dict{Index{d,Vector{N}},N}()
  return Settings{d,V,N,T,I,G,F,typeof(generatorState)}(m0,maxK,gamma,splitting,Z,indexset,numberGenerator,
    gaussianFieldSampler,sampleFunction,showInfo,useTime,generatorState)
end

# reset settings
function reset{S<:Settings}(settings::S)
  reset(settings.numberGenerator)
end

# sampler type
type Sampler{d,N<:Integer,S<:Settings,D<:Dict}
  K::N # keeps track of the size of the index set
  settings::S # simulation details
  samples::D # stores all the samples
end

# outer constructor
function createSampler{S}(settings::S)
  d = ndims(settings)
  N = inttype(settings)
  samples = Dict{Index{d,Vector{N}},Array{Float64,3}}()
  return Sampler{d,Int64,S,typeof(samples)}(0, settings, samples)
end

# convenience method for sampler kind
inttype{d,N,S,D}(sampler::Sampler{d,S,D,N}) = N

ndims{d}(sampler::Sampler{d}) = d

kind(sampler::Sampler) = sampler.settings.indexset.kind

# reset sampler
function reset{S<:Sampler}(sampler::S)
  reset(sampler.settings)
  d = ndims(sampler)
  N = inttype(sampler)
  sampler.samples = Dict{Index{d,Vector{N}},Array{Float64,3}}()
  sampler.K = 0
end

# sample nbOfSamples from the Sampler at the given index
function sample{S<:Sampler,N<:Integer,I<:Index}(sampler::S, nbOfSamples::N, index::I)

  nshift = nshifts(sampler.settings.numberGenerator)
  Z = sampler.settings.Z
  T = typeof(sampler.settings.gamma)

  # print some info to screen
  if isa(sampler.settings.numberGenerator, MCgenerator)
    if ndims(index) == 1
      !sampler.settings.showInfo || println(">> Taking $(nbOfSamples) samples at level $(index[1])...")
    else
      !sampler.settings.showInfo || println(">> Taking $(nbOfSamples) samples at $(index)...")
    end
  else
    if ndims(index) == 1
      !sampler.settings.showInfo || println(">> Taking $(nshift)x$(nbOfSamples) samples at level $(index[1])...")
    else
      !sampler.settings.showInfo || println(">> Taking $(nshift)x$(nbOfSamples) samples at $(index)...")
    end
  end

  # ask state of generator
  state = sampler.settings.generatorState
  if !haskey(state,index)
    state[index] = 0
  end
  s = state[index]

  # take the samples
  samples::Array{T,3} = @parallel (vcat) for i = s+1:s+nbOfSamples

    # generate "random" numbers
    XI = getPoint(sampler.settings.numberGenerator,i)

    # solve
    mySample = zeros(T,1,nshift,Z)::Array{T,3}
    for q = 1:nshift
      xi = XI[:,q]::Vector{T}
      Qtot = 0. # local sum
      j = copy(index)::I
      p = 1
      while j != -1
        Q = sampler.settings.sampleFunction(xi,j,sampler.settings)
        if mod(diff(j,index),2) == 0
          Qtot += Q
        else
          Qtot -= Q
        end
        j = difference(p,j,index)
      end
      mySample[1,q,:] = Qtot
    end
    mySample
  end

  # set state of generator
  state[index] = s+nbOfSamples

  # add samples to the sampler
  if isa(sampler.settings.numberGenerator, MCgenerator)
    samples = permutedims(samples,[2,1,3]) # permute dims if MC
  end

  if !haskey(sampler.samples,index)
    setindex!(sampler.samples,samples,index)
  else
    sampler.samples[index] = isa(sampler.settings.numberGenerator, MCgenerator) ?
      hcat(sampler.samples[index],samples) : vcat(sampler.samples[index],samples)
  end

  return Void

end
