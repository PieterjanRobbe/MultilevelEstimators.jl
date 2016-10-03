# settings type
type Settings{d,V<:AbstractVector,W<:AbstractVector,N<:Integer,T<:AbstractFloat,I<:IndexSet,G<:NumberGenerator,F<:GaussianFieldSampler,D<:Dict}
  indexset::I # type of index set
  numberGenerator::G # type of point generator
  sampleFunction::Function # quantity of interest
  m0::V # coarsest mesh size
  maxL::N # max index number
  γ::W # solver complexity
  splitting::T # desired MSE splitting
  Z::N # number of QOI's to determine
  Nstar::N # number of warm-up samples
  gaussianFieldSampler::F # type of Gaussian field sampler
  showInfo::Bool # boolean to decide if something must be printed to screen
  useTime::Bool # boolean to indicate if true simulation time must be used as cost measure
  safety::Bool # for MLMC, if true then it is guaranteed that variance of estimator is less than TOL^2/2
  generatorState::D # save states of number generator
  procMap::D # maps index to number of processors to be used
end

# utilities
inttype{d,V,W,N}(settings::Settings{d,V,W,N}) = N
floattype{d,V,W,N,T}(settings::Settings{d,V,W,N,T}) = T
ndims{d}(settings::Settings{d}) = d

# outer constructor
function createSettings{I,G,V,W,N,T,F}(
  indexset::I,
  numberGenerator::G,
  sampleFunction::Function;
  m0::V = 4*ones(Int64,1),
  maxL::N = 6,
  γ::W = 1.2*ones(Float64,1),
  splitting::T = 0.5,
  Z::N = 1,
  Nstar::N = 16,
  gaussianFieldSampler::F = EmptySampler(),
  showInfo::Bool = true,
  useTime::Bool = false,
  safety = true,
  procMap = Dict{Index,Integer}() )

  # assertions
  @assert isValid(indexset) && isValid(numberGenerator) && all(m0.>0) && maxL >= 0 && 
    all(γ.>0) && splitting >= 0 && splitting <= 1 && Z >= 1 && isValid(gaussianFieldSampler) && Nstar > 0
  
  # get dimenion
  d = ndims(indexset)

  # check coarsest mesh sizes
  length(m0) == d || warn("length of vector m0 with coarsest mesh sizes is not equal to d, I will try to proceed by adapting length(m0).")
  if length(m0) > d
    m0 = m0[1:d]
  elseif length(m0) < d
    while length(m0) < d
      push!(m0,m0[1])
    end
  end

  # check solver complexities
  length(γ) == d || warn("length of vector γ with solver complexities is not equal to d, I will try to proceed by adapting length(γ).")
  if length(γ) > d
    γ = γ[1:d]
  elseif length(γ) < d
    while length(γ) < d
      push!(γ,γ[1])
    end
  end 

  # assert m0[1] is equal to the one set by the gaussianFieldSampler
  if typeof(gaussianFieldSampler) != EmptySampler
    m0[1] == size(gaussianFieldSampler.eigenfunc[1],2) || warn("coarsest mesh size is not equal for gaussianFieldSampler and this Settings, I will try to proceed anyway.")
  end

  # if we use a QMC generator, safety must be turned on
  G <: QMCgenerator ? safety == true : error("when using a QMC generator, safety must be turned on!")

  # dict for point generator states
  generatorState = Dict{Index{d,Vector{N}},N}()

  # compose worker map, if not empty
  if isempty(procMap)
    id_set = ( d == 1 ? createIndexSet(ML,d) : createIndexSet(FT,d) )
    id_x = getIndexSet(id_set,maxL)
    [setindex!(procMap,nprocs(),Index(i)) for i in id_x]
  else
    minimum(values(procMap)) > 0 || error("wrong procMap provided, nprocs must be positive!")
    maximum(values(procMap)) <= nprocs() || error("I only have $(nprocs()) processor(s) available. Lower the requested number of processors in the procMap or start Julia with more processors.")
  end

  return Settings{d,V,W,N,T,I,G,F,typeof(generatorState)}(indexset,numberGenerator,sampleFunction,m0,maxL,γ,splitting,Z,
    Nstar,gaussianFieldSampler,showInfo,useTime,safety,generatorState,procMap)
end

# methods
function isValid(settings::Settings)
  d = ndims(settings)
  N = inttype(settings)
  returnvalue = isValid(settings.indexset) && isValid(settings.numberGenerator) && all(settings.m0.>0) && settings.maxL >= 0 && 
    all(settings.γ.>0) && settings.splitting >= 0 && settings.splitting <= 1 && settings.Z >= 1 && isValid(settings.gaussianFieldSampler) &&
    length(settings.m0) == d && ( typeof(settings.gaussianFieldSampler) == EmptySampler || settings.m0[1] == size(settings.gaussianFieldSampler.eigenfunc[1],2)) &&
    minimum(values(settings.procMap)) > 0 && maximum(values(settings.procMap)) <= nprocs() && typeof(settings.generatorState) == Dict{Index{d,Vector{N}},N} &&
    typeof(settings.procMap) == Dict{Index{d,Vector{N}},N} && length(settings.γ) == d && settings.Nstar  > 0 && (typeof(settings.numberGenerator) <: QMCgenerator ? settings.safety == true : false)
end

function show(io::IO, settings::Settings)
  itype = ndims(settings.indexset) == 1 ? "level" : "index"
  str  = "  ****************************************************** \n"
  str *= "  *                    Settings with                   * \n"
  str *= "  ****************************************************** \n"
  str *= "  \n"
  str *= "  indexset = $(settings.indexset) \n"
  str *= "  numberGenerator = $(settings.numberGenerator) \n"
  str *= "  sampleFunction = $(settings.sampleFunction) \n"
  str *= "  m0 = $(settings.m0) # coarsest mesh sizes (vector of length d) \n"
  str *= "  maxL = $(settings.maxL) # maximum level parameter for index sets \n"
  str *= "  γ= $(settings.γ) # solver complexities (vector of length d) \n"
  str *= "  splitting = $(settings.splitting) # desired splitting between bias and variance (0.5) \n"
  str *= "  Z = $(settings.Z) # number of qoi's \n"
  str *= "  Nstar = $(settings.Nstar) # number of warm-up samples"
  str *= "  gaussianFieldSampler = $(settings.gaussianFieldSampler) \n"
  str *= "  showInfo = $(settings.showInfo) \n"
  str *= "  useTime = $(settings.useTime) # use true simulation time as cost measure \n"
  str *= "  safety = $(settings.safety) # guarantee that variance of estimator < 0.5*TOL^2 \n"
  str *= "  generator states: \n"
  str *= "  -------------------------------------- \n"
  str *= "    $itype             generator state   \n"
  str *= "  -------------------------------------- \n"
  for index in sort(collect(keys(settings.generatorState)))
    str *= "    $(index.indices)                         "[1:18]
    str *= "    $(settings.generatorState[index])\n"
  end
  str *= "  processor map: \n"
  str *= "  -------------------------------------- \n"
  str *= "    $itype             processor map     \n"
  str *= "  -------------------------------------- \n"
  for index in sort(collect(keys(settings.procMap)))
    str *= "    $(index.indices)                        "[1:18]
    str *= "    $(settings.procMap[index])\n"
  end
  print(io, str)
end

# reset settings
function reset{S<:Settings}(settings::S)
  d = ndims(settings)
  N = inttype(settings)
  reset(settings.numberGenerator)
  settings.generatorState = Dict{Index{d,Vector{N}},N}()
end

# sampler type
type Sampler{d,N<:Integer,S<:Settings,D<:Dict,E<:Dict}
  L::N # keeps track of the size of the index set
  settings::S # simulation settings
  samples::D # stores all the samples
  times::E # store the times for taking these samples
end

# outer constructor
function createSampler{d,V,W,N,T}(settings::Settings{d,V,W,N,T})
  isValid(settings) || error("incorrect settings specified")
  samples = Dict{Index{d,Vector{N}},Array{T,3}}()
  times = Dict{Index{d,Vector{N}},T}()
  return Sampler{d,N,typeof(settings),typeof(samples),typeof(times)}(0, settings, samples, times)
end

# utilities
inttype(sampler::Sampler) = inttype(sampler.settings)
floattype(sampler::Sampler) = floattype(sampler.settings)
ndims{d}(sampler::Sampler{d}) = d
kind(sampler::Sampler) = kind(sampler.settings.indexset)

# methods
function isValid(S::Sampler)
  d = ndims(sampler.settings)
  N = inttype(sampler.settings)
  T = floattype(sampler.settings)
  returnvalue = L > 0 && isValid(sampler.settings) && typeof(sampler.samples) == Dict{Index{d,Vector{N}},Array{T,3}} &&
  typeof(sampler.times) == Dict{Index{d,Vector{N}},T}
end

function show(io::IO, sampler::Sampler)
  itype = ndims(sampler.settings.indexset) == 1 ? "level" : "index"
  str  = "  ****************************************************** \n"
  str *= "  *                    Sampler with                    * \n"
  str *= "  ****************************************************** \n"
  str *= "  \n"
  str *= "  L = $(sampler.L) # current level parameter \n"
  substr = isValid(sampler.settings) ? "*valid*" : "*invalid* !!!"
  str *= "  settings = $(substr) \n"
  str *= "  samples: \n"
  str *= "  ------------------------------------------------------ \n"
  str *= "    $itype           sample size           time/sample   \n"
  str *= "  ------------------------------------------------------ \n"
  for index in sort(collect(keys(sampler.samples)))
    str *= "    $(index.indices)                         "[1:16]
    str *= "    $(size(sampler.samples[index]))                   "[1:22]
    str *= @sprintf("    %0.6e\n",sampler.times[index])
  end
  print(io, str)
end

function reset{S<:Sampler}(sampler::S)
  reset(sampler.settings)
  d = ndims(sampler)
  N = inttype(sampler)
  T = floattype(sampler)
  sampler.samples = Dict{Index{d,Vector{N}},Array{T,3}}()
  sampler.times = Dict{Index{d,Vector{N}},T}()
  sampler.L = 0
end

#
# Main methods for sampling from the sampler
#

# sample nbOfSamples from the Sampler at the given index
function sample{S<:Sampler,N<:Integer,I<:Index}(sampler::S, nbOfSamples::N, index::I)

  # print some info to screen
  itype = ndims(sampler.settings.indexset) == 1 ? "level" : "index"
  if isa(sampler.settings.numberGenerator, MCgenerator)
      !sampler.settings.showInfo || println(">> taking $(nbOfSamples) samples at "*itype*" $(index.indices)...")
  else
      !sampler.settings.showInfo || println(">> taking $(nshifts(sampler.settings.numberGenerator))x$(nbOfSamples) samples at "*itype*" $(index.indices)...")
  end

  # ask state of generator
  state = sampler.settings.generatorState
  if !haskey(state,index)
    state[index] = 0
  end
  gen_state = state[index]

  # total number of workers to be used
  p = max(1,get(sampler.settings.procMap,index,Void) - 1) # nworkers() = nprocs() - 1

  # take the samples
  r = RemoteRef{Channel{Any}}[]
  for pid = workers()[1]:workers()[1]+p-1
    push!(r,remotecall(pid,()->sample_p(sampler,nbOfSamples,index,gen_state)))
  end

  samples = reduce(vcat,map(fetch,r))
  @assert size(samples,1) == nbOfSamples

  # set state of generator
  state[index] += nbOfSamples

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

# defines subfunction for pmap operation
function sample_p{S<:Sampler,N<:Integer,I<:Index}(sampler::S, n::N, index::I, gen_state::N)

  # aliases
  nshift = nshifts(sampler.settings.numberGenerator)
  Z = sampler.settings.Z
  T = floattype(sampler)

  # block distribution
  p = max(1,get(sampler.settings.procMap,index,Void) - 1) # nworkers() = nprocs() - 1
  s = max(1,myid() - 1) # processor id, range(2,p+1)
  b = convert(N,ceil(n/p)) # block size

  # preallocate samples
  samples = zeros(T,max(0,min(n,s*b) - (s - 1)*b),nshift,Z)

  for i = gen_state + (s - 1)*b + 1 : gen_state + min(n,s*b)

    # generate "random" numbers
    XI = getPoint(sampler.settings.numberGenerator,i)

    # solve
    mySample = zeros(T,1,nshift,Z)::Array{T,3}
    for q = 1:nshift
      xi = XI[:,q]::Vector{T}
      Qtot = 0. # local sum
      j = copy(index)::I
      while j != -1
        Q = sampler.settings.sampleFunction(xi,copy(j),sampler.settings) # ! copy
        if mod(diff(j,index),2) == 0
          Qtot += Q
        else
          Qtot -= Q
        end
        j = difference(j,index)
      end
      mySample[1,q,:] = Qtot
    end
    samples[i - gen_state - (s - 1)*b,:,:] = mySample
  end

  return samples
end
