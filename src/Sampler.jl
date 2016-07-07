# simulation settings
type Settings{T,N} # only mutable type is passed by reference!
  m0::Vector{N} # coarsest mesh size
  maxK::N # max index number
  gamma::T # solver complexity
  splitting::T # desired MSE splitting
  Z::N # number of QOI's to determine
  indexset::IndexSet # type of index set
  numberGenerator::NumberGenerator # type of point generator
  gaussianFieldSampler::GaussianFieldSampler # type of Gaussian field sampler
  sampleFunction::Function # quantity of interest
  showInfo::Bool # boolean to decide if something must be printed to screen
  useTime::Bool # boolean to indicate if true simulation time must be used as cost measure
end

# convenience constructor for settings
function Settings{d,T,N}(indexset::IndexSet{d}, numberGenerator::NumberGenerator, gaussianFieldSampler::GaussianFieldSampler, sampleFunction::Function;
  m0::Vector{N} = 4*ones(Int64,d), maxK::N = 6, gamma::T = 1.4, splitting::T = 0.5, Z::N = 1, showInfo::Bool = true, useTime::Bool = false)
  return Settings(m0,maxK,gamma,splitting,Z,indexset,numberGenerator,gaussianFieldSampler,sampleFunction,showInfo,useTime)
end

function reset!{T,N}(settings::Settings{T,N})
  reset!(settings.numberGenerator)
end

# sampler type
type Sampler{d,T<:AbstractFloat,N<:Integer}
  settings::Settings{T,N} # simulation details
  samples::Dict{Index,Array{T,3}} # stores all the samples
  K :: N # keeps track of the size of the index set
end

# convenience method for sampler kind
kind{d,T,N}(sampler::Sampler{d,T,N}) = sampler.settings.indexset.kind

function reset!{T,N}(sampler::Sampler{T,N})
  reset!(sampler.settings)
  sampler.samples = Dict{Index,Array{T,1}}()
  sampler.K = 0
end

# constructor
function createSampler{T,N}(d::N,settings::Settings{T,N})
  return Sampler{d,T,N}(settings, Dict{Index,Array{T,1}}(), 0, )
end

# sample nbOfSamples from the Sampler at the given index
function sample{d,T,N}(sampler::Sampler{d,T,N}, nbOfSamples::N, index::Index{d})

  nshift = nshifts(sampler.settings.numberGenerator)
  Z = sampler.settings.Z

  # print some info to screen
  if isa(sampler.settings.numberGenerator, MCgenerator)
    !sampler.settings.showInfo || println(">> Taking $(nbOfSamples) samples at $(index)...")
  else
    !sampler.settings.showInfo || println(">> Taking $(nshift)x$(nbOfSamples) samples at $(index)...")
  end

  # take the samples
  samples = @parallel (vcat) for i = 1:nbOfSamples

    # generate "random" numbers
    XI = next(sampler.settings.numberGenerator,index)

    # solve
    mySample = zeros(T,1,nshift,Z)
    for q = 1:nshift
      xi = XI[:,q]
      Qtot = 0 # local sum
      j = copy(index)
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

end
