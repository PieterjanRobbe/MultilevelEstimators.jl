# sampler type
type Sampler{d,I<:IndexSet,G<:NumberGenerator,F,D1<:Dict,D2<:Dict,D3<:Dict,UType}

  # settings that must be provided by the user
  indexSet::I                   # type of index set
  numberGenerator::G            # type of point generator
  sampleFunction::Function      # quantity of interest (qoi)

  # settings that can be altered by the user
# coarseGrids::Vector{Float64}  # coarsest mesh sizes (vector of length d) (obsolete ???)
  maxL::Int64                   # max index number
  γ::Vector{Float64}            # solver complexity (vector of length d)
# splitting::Float64            # desired MSE splitting (obsolete ???)
  Z::Int64                      # number of qoi's to determine
  Nstar::Int64                  # number of warm-up samples
  gaussianFieldSampler::F       # type of Gaussian field sampler
  useTime::Bool                 # use true simulation time as cost measure
  safety::Bool                  # guarantee that variance of estimator is less than θ^2*TOL^2
  continuate::Bool              # use Bayesian inversion for parameter estimation
  nTOL::Int64                   # number of iterations in the continuation algorithm
  k::Tuple{Float64,Float64}     # trust parameters for normal-gamma prior in continuation algorithm
  showInfo::Bool                # print diagnostics to screen
  ioStream::IO                  # where to write output to (file or std out)
  storeSamples0::Bool           # store non-difference samples for problem analysis
  procMap::D1                   # maps index to number of processors to be used
  userType::UType

  # unchangeable
  generatorStates::D1           # save states of number generator
  samples::D2                   # stores all the difference samples
  samples0::D2                  # stores all the non-difference samples
  times::D3                     # store the time / index
end

# utilities
ndims{d}(sampler::Sampler{d}) = d

# constructor
function setup{S<:AbstractString}(dict::Dict{S,Any})

  settings = copy(dict)

  if !haskey(settings,"indexSet")
    error("I need an indexSet (SL/ML/FT/TD/HC/AD) before I can do anything!")
  else
    indexSet = settings["indexSet"]
    delete!(settings,"indexSet")
    if !(typeof(indexSet) <: IndexSet)
      error("incorrect indexSet specified!")
    end
  end
  d = ndims(indexSet)

  if !haskey(settings,"numberGenerator")
    error("I need a numberGenerator before I can do anything!")
  else
    numberGenerator = settings["numberGenerator"]
    delete!(settings,"numberGenerator")
    if !(typeof(numberGenerator) <: NumberGenerator)
      error("incorrect numberGenerator specified!")
    end
  end

  if !haskey(settings,"sampleFunction")
    error("I need a sampleFunction (aka qoi) before I can do anything!")
  else
    sampleFunction = settings["sampleFunction"]
    delete!(settings,"sampleFunction")
    if !(typeof(sampleFunction) <: Function)
      error("incorrect sampleFunction specified!")
    end
  end

  # if haskey(settings,"coarseGrids")
  #   if length(coarseGrids) > d
  #     warn("length of vector coarseGrids with coarsest mesh sizes is larger than d, but I will try to proceed by adapting length(coarseGrids).")
  #     coarseGrids = coarseGrids[1:d]
  #   elseif length(coarseGrids) < d
  #     warn("length of vector coarseGrids with coarsest mesh sizes is smaller than d, but I will try to proceed by adapting length(coarseGrids).")
  #     while length(coarseGrids) < d
  #       push!(coarseGrids,coarseGrids[1])
  #     end
  #   end
  # else
  #   coarseGrids = 4*ones(Int64,d)
  # end

  if haskey(settings,"maxL")
    maxL = settings["maxL"]
    delete!(settings,"maxL")
    if !(typeof(maxL) == Int64)
      error("incorrect maxL specified!")
    elseif maxL < 0
      error("maximum level parameter maxL must be positive!")
    end
  else
    maxL = 6
  end
  if ( typeof(indexSet) == SL && maxL != 0 )
    warn("for SL simulation, maxL must equal zero. I will try to proceed with maxL=0")
    maxL = 0
  end

  if haskey(settings,"γ")
    γ = settings["γ"]
    # delete!(settings,"γ")
    if !(typeof(γ) == Vector{Float64})
      error("incorrect solver complexities γ specified!")
    elseif length(γ) > d
      warn("length of vector γ with solver complexities is larger than d, but I will try to proceed by adapting length(γ).")
      γ = γ[1:d]
    elseif length(γ) < d
      warn("length of vector γ with solver complexities is smaller than d, but I will try to proceed by adapting length(γ).")
      while length(γ) < d
        push!(γ,γ[1])
      end
    elseif !all(γ.>0)
      error("solver complexities γ must be positive!")
    elseif any(γ.<1)
      warn("one of the solver complexities γ < 1, are you sure this is correct?")
    end
  else
    γ = 1.5*ones(Float64,d)
  end

  # if haskey(settings,"splitting")
  #   splitting = settings["splitting"]
  #   delete!(settings,"splitting")
  #   if !(typeof(splitting) == Float64)
  #     error("incorrect splitting specified!")
  #   elseif (splitting < 0 || splitting > 1)
  #     error("splitting parameter must be between 0 and 1")
  #   end
  # else
  #   splitting = 0.5
  # end

  if haskey(settings,"Z")
    Z = settings["Z"]
    delete!(settings,"Z")
    if !(typeof(Z) == Int64)
      error("incorrect number of qoi's Z specified!")
    elseif Z < 1
      error("number of qoi's Z must be positive!")
    end
  else
    Z = 1
  end

  if haskey(settings,"Nstar")
    Nstar = settings["Nstar"]
    delete!(settings,"Nstar")
    if !(typeof(Nstar) == Int64)
      error("incorrect number of warm-up samples Nstar specified!")
    elseif Nstar < 1
      error("number of warm-up samples Nstar must be positive!")
    end
  else
    Nstar = convert(Int64,ceil(32/nshifts(numberGenerator)))
  end

  if haskey(settings,"gaussianFieldSampler")
    gaussianFieldSampler = settings["gaussianFieldSampler"]
    delete!(settings,"gaussianFieldSampler")
    if !(typeof(gaussianFieldSampler) <: GaussianFieldSampler || (typeof(gaussianFieldSampler) <: Array && eltype(gaussianFieldSampler) <: GaussianFieldSampler ) )
      error("incorrect Gaussian field sampler gaussianFieldSampler specified!")
    elseif typeof(gaussianFieldSampler) == KLExpansion && maxL+1 < length(gaussianFieldSampler.eigenfunc)
      warn("you asked for $(maxL+1) levels, but the KL expansion only provides $(length(gaussianFieldSampler.eigenfunc)) eigenfunctions. I will try to continue anyway.")
    elseif typeof(gaussianFieldSampler) <: Array
      for i in 1:length(gaussianFieldSampler)
        if typeof(gaussianFieldSampler[i]) == KLExpansion && maxL+1 < length(gaussianFieldSampler[i].eigenfunc)
          warn("you asked for $(maxL+1) levels, but the $(i)-th KL expansion only provides $(length(gaussianFieldSampler.eigenfunc)) eigenfunctions. I will try to continue anyway.")
        end
      end
    end
    # if coarseGrids[1] != size(gaussianFieldSampler.eigenfunc[1],2)
    #   warn("coarsest mesh size in coarseGrids is not equal to the mesh size of the eigenfunctions in gaussianFieldSampler, but I will try to proceed anyway.")
    # end
  else
    gaussianFieldSampler = EmptySampler()
  end

  if haskey(settings,"useTime")
    useTime = settings["useTime"]
    delete!(settings,"useTime")
    if !(typeof(useTime) == Bool)
      error("incorrect boolean useTime specified!")
    end
  else
    useTime = false
  end
  # check that either useTime is false and γ was provided or useTime is true and a cost model γ is provided
  if useTime == false && !haskey(settings,"γ")
    warn("useTime is set to false, but no cost model γ was provided. I will use the standard cost model γ= 1.5*ones(d) and try to continue.")
  elseif useTime == true && haskey(settings,"γ")
    warn("useTime is set to true, but also a cost model γ was provided. I will use the true simulation times and try to continue.")
  end
  delete!(settings,"γ")

  if haskey(settings,"safety")
    safety = settings["safety"]
    delete!(settings,"safety")
    if !(typeof(safety) == Bool)
      error("incorrect boolean safety specified!")
    end
  else
    safety = true
  end
  if typeof(numberGenerator) <: QMCgenerator && safety == false
    warn("when using a QMC generator, safety must be turned on!")
    safety = true
  end

  if haskey(settings,"continuate")
    continuate = settings["continuate"]
    delete!(settings,"continuate")
    if !(typeof(continuate) == Bool)
      error("incorrect boolean continuate specified!")
    end
  else
    continuate = false
  end

  if haskey(settings,"nTOL")
    nTOL = settings["nTOL"]
    delete!(settings,"nTOL")
    if !(typeof(nTOL) == Int64)
      error("incorrect number of continuation tolerances nTOL specified!")
    elseif nTOL < 0
      error("number of continuation tolerances nTOL must be positive!")
    elseif !continuate
      warn("number of continuation tolerances nTOL was provided, but continuation is set to false. I will ignore entry 'nTOL'")
    end
  else
    nTOL = 10
  end

  if haskey(settings,"k")
    k = settings["k"]
    delete!(settings,"k")
    if !(typeof(k) == Tuple{Float64,Float64})
      error("incorrect trust parameters k specified!")
    elseif length(k) != 2
      error("I can only accept 2 trust parameters k0 and k1")
    elseif k[1]<0 || k[2]<0
      error("trust parameters k must be positive")
    elseif !continuate
      warn("trust parameters k where provided, but continuation is set to false. I will ignore entry 'k'")
    end
  else
    k = (0.1,0.1)
  end

  if haskey(settings,"showInfo")
    showInfo = settings["showInfo"]
    delete!(settings,"showInfo")
    if !(typeof(showInfo) == Bool)
      error("incorrect boolean showInfo specified!")
    end
  else
    showInfo = true
  end

  if haskey(settings,"ioStream")
    ioStream = settings["ioStream"]
    if !(typeof(ioStream) == IO)
      error("incorrect ioStream specified!")
    end
  else
    ioStream = STDOUT
  end

  if haskey(settings,"storeSamples0")
    storeSamples0 = settings["storeSamples0"]
    delete!(settings,"storeSamples0")
    if !(typeof(storeSamples0) == Bool)
      error("incorrect boolean storeSamples0 specified!")
    end
  else
    storeSamples0 = false
  end

  if haskey(settings,"procMap")
    procMap = settings["procMap"]
    delete!(settings,"procMap")
    if !(typeof(ioStream) == Dict{Index{d,Vector{Int64}},Int64})
      error("incorrect procMap specified!")
    elseif minimum(values(procMap)) <= 0
      error("wrong procMap provided, nprocs must be positive!")
    elseif maximum(values(procMap)) > nprocs()
      error("I only have $(nprocs()) processor(s) available. Lower the requested number of processors in the procMap or start Julia with more processors.")
    end
  else
    procMap = Dict{Index{d,Vector{Int64}},Int64}()
    dummyIndexSet = typeof(indexSet) <: AD ? TD(d) : indexSet
    id_x = getIndexSet(dummyIndexSet,min(20,maxL)) # hack
    [setindex!(procMap,nprocs(),i) for i in id_x]
  end

  if haskey(settings,"userType")
    userType = settings["userType"]
    delete!(settings,"userType")
  else
    userType = Void()
  end

  if !isempty(settings)
    str = ""
    for s in keys(settings)
      str *= s*", "
    end
    str = str[1:end-2]
    error("could not parse the following entries: "*str)
  end

  generatorStates = Dict{Index{d,Vector{Int64}},Int64}()

  samples = Dict{Index{d,Vector{Int64}},Array{Float64,3}}()

  samples0 = Dict{Index{d,Vector{Int64}},Array{Float64,3}}()

  times = Dict{Index{d,Vector{Int64}},Float64}()
  
	return Sampler{d,typeof(indexSet),typeof(numberGenerator),typeof(gaussianFieldSampler),typeof(generatorStates),typeof(samples),typeof(times),typeof(userType)}(
    indexSet, numberGenerator, sampleFunction, maxL, γ, Z, Nstar, gaussianFieldSampler, useTime,
      safety, continuate, nTOL, k, showInfo, ioStream, storeSamples0, procMap, userType, generatorStates, samples, samples0, times)
end

function show(io::IO, sampler::Sampler)
  itype = ndims(sampler) == 1 ? "level" : "index"
  str  = "********************************************************************************\n"
  str *= "*                                SAMPLER STATUS                                *\n"
  str *= "********************************************************************************\n"
  str *= "indexSet        | $(sampler.indexSet) \n"
  str *= "numberGenerator | $(sampler.numberGenerator) \n"
  str *= "sampleFunction  | $(sampler.sampleFunction) \n"
# str *= "coarseGrids     | $(sampler.coarseGrids) # coarsest mesh sizes (vector of length d) \n"
  str *= "maxL            | $(sampler.maxL) # maximum level parameter for index sets \n"
  str *= "γ               | $(sampler.γ) # solver complexities (vector of length d) \n"
# str *= "splitting       | $(sampler.splitting) # desired splitting between bias and variance (0.5) \n"
  str *= "Z               | $(sampler.Z) # number of qoi's \n"
  str *= "Nstar           | $(sampler.Nstar) # number of warm-up samples \n"
  str *= "gaussianFiel... | $(sampler.gaussianFieldSampler) \n"
  str *= "useTime         | $(sampler.useTime) # use true simulation time as cost measure \n"
  str *= "safety          | $(sampler.safety) # guarantee that variance of estimator < θ*TOL^2 \n"
  str *= "continuate      | $(sampler.continuate) # do Bayesian parameter estimation \n"
  str *= "k               | $(sampler.k) # trust parameters in normal-gaussian prior \n"
  str *= "showInfo        | $(sampler.showInfo) \n"
  str *= "storeSamples0   | $(sampler.storeSamples0) \n"
  str *= "ioStream        | $(sampler.ioStream) \n"
  str *= "userType        | $(sampler.userType) \n"
  str *= "--> processor map: \n"
  str *= "-------------------------------------- \n"
  str *= "  $itype             processor map     \n"
  str *= "-------------------------------------- \n"
  for index in sort(collect(keys(sampler.procMap)))
    str *= "  $(index.indices)                        "[1:18]
    str *= "  $(sampler.procMap[index])\n"
  end
  str *= "--> generator states: \n"
  str *= "-------------------------------------- \n"
  str *= "  $itype             generator state   \n"
  str *= "-------------------------------------- \n"
  for index in sort(collect(keys(sampler.generatorStates)))
    str *= "  $(index.indices)                         "[1:18]
    str *= "  $(sampler.generatorStates[index])\n"
  end
  dir = isa(sampler.numberGenerator,MCgenerator) ? 2 : 1
  str *= "--> samples: \n"
  str *= "------------------------------------------------------ \n"
  str *= "  $itype           sample size           time/sample   \n"
  str *= "------------------------------------------------------ \n"
  for index in sort(collect(keys(sampler.samples)))
    nsamples = size(sampler.samples[index],dir)
    str *= "  $(index.indices)                         "[1:16]
    str *= "  $(size(sampler.samples[index]))                   "[1:22]
    str *= @sprintf("  %0.6e\n",sampler.times[index]/nsamples)
  end
  print(io, str)
end

# save current state of the sampler for later reinitialisation
function save(sampler::Sampler)
  # make directory structure
  current_dir = pwd()
  dir_name = "data_"*Dates.format(now(),"dd_mm_yy_HH_MM")
  mkdir(dir_name)
  cd(dir_name)

  # print sampler info to info.txt
  file_handle = open("info.txt","a")
  show(file_handle, sampler)
  close(file_handle)

  # save all samples
  dir_name = "samples"
  mkdir(dir_name)
  cd(dir_name)
  for idx in keys(sampler.samples)
    dir_name = @sprintf("%s",idx.indices)
    mkdir(dir_name)
    cd(dir_name)
    writedlm("gen_state.txt",sampler.generatorStates[idx])
    writedlm("time.txt",sampler.times[idx])
    for z in 1:sampler.Z
      writedlm("Z$(z).txt",sampler.samples[idx][:,:,z])
    end
    cd("..")
  end
  cd("..")

  # save all samples0 (non-difference samples)
  if sampler.storeSamples0
    dir_name = "samples0"
    mkdir(dir_name)
    cd(dir_name)
    for idx in keys(sampler.samples)
      dir_name = @sprintf("%s",idx.indices)
      mkdir(dir_name)
      cd(dir_name)
      for z in 1:sampler.Z
        writedlm("Z$(z).txt",sampler.samples0[idx][:,:,z])
      end
      cd("..")
    end
    cd("..")
  end

  # print mimc output table to current_state.txt
  d = ndims(sampler)
  E     = Dict{Index{d,Vector{Int64}},Vector{Float64}}()
  V     = Dict{Index{d,Vector{Int64}},Vector{Float64}}()
  W     = Dict{Index{d,Vector{Int64}},Float64}() # dict for costs
  for index in keys(sampler.samples)
    E[index] = squeeze(mean(sampler.samples[index],(1,2)),(1,2))
    V[index] = squeeze(var(mean(sampler.samples[index],1),2),(1,2))
    W[index] =  sampler.useTime ? sampler.times[index] : prod(2.^(index.indices.*sampler.γ))
  end
  str  = "-------------------------------------------------------------------------------- \n"
  itype = ndims(sampler.indexSet) == 1 ? "level" : "index"
  str *= "  "*itype*"       E              V               N               W               \n"
  str *= "-------------------------------------------------------------------------------- \n"
  for index in sort(Set(collect(keys(E))))
    str *= ("  $(index.indices)            "[1:13])
    str *= @sprintf("%12.5e",maximum(E[index]))
    str *= @sprintf("    %0.6e",maximum(V[index]))
    str *= @sprintf("    %d               ",prod(size(sampler.samples[index])[1:2]))[1:16]
    str *= @sprintf("    %0.6e \n",W[index])
  end
  file_handle = open("current_state.txt","a")
  @printf(file_handle, "%s", str)
  close(file_handle)

  cd("..")
end

# reinitialize sampler
function load(sampler::Sampler,dir::AbstractString)
  isdir(dir) || error("$(dir) is not a directory")
  cd(dir*"/samples")
  for idx_dir in filter!(r",", readdir())
    cd(idx_dir)
    if !(idx_dir[1] == ".") # ignore hidden files/folders
      idx = Index(idx_dir)
      setindex!(sampler.generatorStates, readdlm("gen_state.txt")[1], idx)
      setindex!(sampler.times, readdlm("time.txt")[1], idx)
      # read Z1 to infer size of samples
      (n,m) = size(readdlm("Z1.txt"))
      l = length(filter!(r"^Z", readdir()))
      samples = zeros(Float64,n,m,l)
      for z in 1:l
        samples[:,:,z] = readdlm("Z$(z).txt")
      end
      setindex!(sampler.samples, samples, idx)
    end
    cd("..")
  end
  cd("../..")
  if sampler.storeSamples0
    cd(dir*"/samples0")
    for idx_dir in filter!(r",", readdir())
      cd(idx_dir)
      if !(idx_dir[1] == ".") # ignore hidden files/folders
        idx = Index(idx_dir)
        # read Z1 to infer size of samples0
        (n,m) = size(readdlm("Z1.txt"))
        l = length(filter!(r"^Z", readdir()))
        samples0 = zeros(Float64,n,m,l)
        for z in 1:l
          samples0[:,:,z] = readdlm("Z$(z).txt")
        end
        setindex!(sampler.samples0, samples0, idx)
      end
      cd("..")
    end
  end
  cd("../..")
end

# reset sampler
function reset(sampler::Sampler)
  d = ndims(sampler)
  sampler.generatorStates = Dict{Index{d,Vector{Int64}},Int64}()
  sampler.L = 0
  sampler.samples = Dict{Index{d,Vector{Int64}},Array{Float64,3}}()
  sampler.samples0 = Dict{Index{d,Vector{Int64}},Array{Float64,3}}()
  sampler.times = Dict{Index{d,Vector{Int64}},Float64}()
end

#
# Main methods for sampling from the sampler
#

# sample nbOfSamples from the Sampler at the given index
function sample{N<:Integer,I<:Index}(sampler::Sampler, nbOfSamples::N, index::I)

  # print some info to screen
  itype = ndims(sampler.indexSet) == 1 ? "level" : "index"
  idcs = ndims(sampler) == 1 ? 1 : 1:ndims(sampler)
  if isa(sampler.numberGenerator, MCgenerator)
      !sampler.showInfo || println(">> taking $(nbOfSamples) samples at "*itype*" $(index.indices[idcs])...")
  else
      !sampler.showInfo || println(">> taking $(nshifts(sampler.numberGenerator))x$(nbOfSamples) samples at "*itype*" $(index.indices[idcs])...")
  end

  # ask state of generator
  state = sampler.generatorStates
  if !haskey(state,index)
    state[index] = 0
  end
  #gen_state = state[index]
	@debug println("current state of point generator is $(state[index])")
  
	# total number of workers to be used
  p = max(1,get(sampler.procMap,index,Void) - 1) # nworkers() = nprocs() - 1
	wp = WorkerPool(collect(1:p)+1)	# workerpool
 	@debug println("running with $(p) workers: $(wp)")

	# function definition
	f(i) = take_a_sample(index,state[index]+i,sampler)

	# parallel execution
	r = pmap(wp, f, 1:nbOfSamples)




#    sample_p(sampler,nbOfSamples,index,gen_state)
  # take the samples
  #r = Future[]
  #myfunc = () -> sample_p(sampler,nbOfSamples,index,gen_state)
  #startiter = workers()[1]
  #for pid = startiter:startiter+p-1
#	  println(pid)
 #   push!(r,remotecall(myfunc,pid))
  #  	println("pushed")
  #end

  #r = [@spawnat pid sample_p(sampler,nbOfSamples,index,gen_state) for pid in startiter:startiter+p-1]


#  println("start...")
  #r = Vector{Array{Float64,4}}(nworkers())
  #@sync begin
  #for (idx, pid) in enumerate(workers())
  #  @async r[idx] = remotecall_fetch(myfunc, pid)
  #end
  #end
#  println("stop...")
  
  #r = Vector{Array{Float64,4}}(1)
	#r[1] = myfunc()

  comb_samples = reduce(vcat,r)::Array{Float64,4}
  @assert size(comb_samples,1) == nbOfSamples

  samples = comb_samples[:,:,:,1]
  samples0 = sampler.storeSamples0 ? comb_samples[:,:,:,2] : zeros(0,0,0)

  # set state of generator
  state[index] += nbOfSamples
	@debug println("current state of point generator is $(state[index])")

  # add samples to the sampler
  if isa(sampler.numberGenerator, MCgenerator)
    samples = permutedims(samples,[2,1,3]) # permute dims if MC
    samples0 = permutedims(samples0,[2,1,3])
  end

  if !haskey(sampler.samples,index)
    setindex!(sampler.samples,samples,index)
    setindex!(sampler.samples0,samples0,index)
  else
    sampler.samples[index] = isa(sampler.numberGenerator, MCgenerator) ?
      hcat(sampler.samples[index],samples) : vcat(sampler.samples[index],samples)
    sampler.samples0[index] = isa(sampler.numberGenerator, MCgenerator) ?
      hcat(sampler.samples0[index],samples0) : vcat(sampler.samples0[index],samples0)
  end
  return Void
end


function take_a_sample{N<:Integer,I<:Index}(index::I, i::N, sampler::Sampler)

	# getting some info from the sampler 
	nshift = nshifts(sampler.numberGenerator)
  Z = sampler.Z
  T = Float64

  # generate "random" numbers
  XI = getPoint(sampler.numberGenerator,i)

  # solve
  mySample = zeros(T,1,nshift,Z,2)::Array{T,4}
  for q = 1:nshift
    xi = XI[:,q]::Vector{T}
    Qtot = 0. # local sum
    Q0 = 0.
    j = copy(index)::I
    while j != -1
      Q = sampler.sampleFunction(xi,copy(j),sampler) # ! copy
      if j == index
        Q0 = Q
      end
      if mod(diff(j,index),2) == 0
        Qtot += Q
      else
        Qtot -= Q
      end
      j = difference(j,index)
    end
    mySample[1,q,:,1] = Qtot
    mySample[1,q,:,2] = Q0
  end
  return mySample
	
	#  println("donecomp")
	return samples
end




# defines subfunction for pmap operation, but now returns samples and samples0 (not differenced)
function sample_p{N<:Integer,I<:Index}(sampler::Sampler, n::N, index::I, gen_state::N)
#	println("begincomp")
	
	nshift = nshifts(sampler.numberGenerator)
  Z = sampler.Z
  T = Float64

  # block distribution
  p = max(1,get(sampler.procMap,index,Void) - 1) # nworkers() = nprocs() - 1
  s = max(1,myid() - 1) # processor id, range(2,p+1)
  b = convert(N,ceil(n/p)) # block size

  # preallocate samples
  samples = zeros(T,max(0,min(n,s*b) - (s - 1)*b),nshift,Z,2)

 # println(max(0,min(n,s*b) - (s - 1)*b))


 for i = gen_state + (s - 1)*b + 1 : gen_state + min(n,s*b)

    # generate "random" numbers
    XI = getPoint(sampler.numberGenerator,i)

    # solve
    mySample = zeros(T,1,nshift,Z,2)::Array{T,4}
    for q = 1:nshift
      xi = XI[:,q]::Vector{T}
      Qtot = 0. # local sum
      Q0 = 0.
      j = copy(index)::I
      while j != -1
        Q = sampler.sampleFunction(xi,copy(j),sampler) # ! copy
        if j == index
          Q0 = Q
        end
        if mod(diff(j,index),2) == 0
          Qtot += Q
        else
          Qtot -= Q
        end
        j = difference(j,index)
      end
      mySample[1,q,:,1] = Qtot
      mySample[1,q,:,2] = Q0
    end
    samples[i - gen_state - (s - 1)*b,:,:,:] = mySample
	end
#  println("donecomp")
return samples
end
