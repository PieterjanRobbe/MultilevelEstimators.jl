# sampler type
type Sampler{d,I<:IndexSet,G<:NumberGenerator,F,D1<:Dict,D2<:Dict,D3<:Dict,D4<:Dict,UType}

  # settings that must be provided by the user
  indexSet::I                   # type of index set
  numberGenerator::G            # type of point generator
  sampleFunction::Function      # quantity of interest (qoi)

  # settings that can be altered by the user
  maxL::Int64                   # max index number
  γ::Vector{Float64}            # solver complexity (vector of length d)
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
  samples::D2                   # all the difference samples
  samples0::D2                  # all the non-difference samples
  T::D3                     		# cumulative time / index
	E::D4													# mean value at each index 
	Vf::D4												# mean of variances at each index
	V::D4													# variance of means at each index (for QMC)
	Vest::D4											# contribution of the variance at each index
	Wst::D3												# standard cost at each index
	W::D3													# true cost at each index used in algorithm, either by runtime or standard cost
	P::D3													# profits at each index
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
	@debug showInfo = true

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
    if !(typeof(procMap) == Dict{Index{d,Vector{Int64}},Int64})
      error("incorrect procMap specified!")
    elseif minimum(values(procMap)) <= 0
      error("wrong procMap provided, nprocs must be positive!")
    elseif maximum(values(procMap)) > nprocs()
      error("I only have $(nprocs()) processor(s) available. Lower the requested number of processors in the procMap or start Julia with more processors.")
    end
  else
    procMap = Dict{Index{d,Vector{Int64}},Int64}()
    dummyIndexSet = typeof(indexSet) <: AD ? TD(d) : indexSet
    id_x = getIndexSet(dummyIndexSet,max(20,maxL)) # hack
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

  T = Dict{Index{d,Vector{Int64}},Float64}()
	E = Dict{Index{d,Vector{Int64}},Vector{Float64}}()
	V = Dict{Index{d,Vector{Int64}},Vector{Float64}}()
	Vf = Dict{Index{d,Vector{Int64}},Vector{Float64}}()
  Wst = Dict{Index{d,Vector{Int64}},Float64}()
  W = Dict{Index{d,Vector{Int64}},Float64}()
	Vest = Dict{Index{d,Vector{Int64}},Vector{Float64}}()
  P = Dict{Index{d,Vector{Int64}},Float64}()
  
	return Sampler{d,typeof(indexSet),typeof(numberGenerator),typeof(gaussianFieldSampler),typeof(generatorStates),typeof(samples),typeof(T),typeof(E),typeof(userType)}(
    indexSet, numberGenerator, sampleFunction, maxL, γ, Z, Nstar, gaussianFieldSampler, useTime,
      safety, continuate, nTOL, k, showInfo, ioStream, storeSamples0, procMap, userType, 
			generatorStates, samples, samples0, T,E,V,Vf,Wst,W,Vest,P)
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

function inspect(sampler)
 	@printf(sampler.ioStream,"%s","--------------------------------------------------------------------------------\n")
 	itype = ndims(sampler.indexSet) == 1 ? "level" : "index"
 	@printf(sampler.ioStream,"%s","  "*itype*"       E              V               N               W              \n")
 	@printf(sampler.ioStream,"%s","--------------------------------------------------------------------------------\n")
 	# TODO store index set in sampler, gives a correct image of all samples for plotting
	for index in sort(Set(keys(sampler.E)))
    str = ("  $(index.indices)            "[1:13])
    str *= @sprintf("%12.5e",maximum(sampler.E[index]))
    str *= @sprintf("    %0.6e",maximum(sampler.Vf[index]))
		str *= @sprintf("    %d               ",prod(size(get(sampler.samples,index,zeros(0,0)))[1:2]))[1:16]
    str *= @sprintf("    %0.6e\n",sampler.W[index])
    @printf(sampler.ioStream,"%s",str)
	end
end

#
# Main methods for sampling from the sampler
#

# compute bias of estimator using the indices in the boundary set
# TODO if we store the current index set in the sampler, we do not need to pass the boundary
function compute_bias{d,N}(sampler, boundary::Set{Index{d,Vector{N}}})
  B = zeros(Float64,sampler.Z)
	for index::Index{d,Vector{N}} in boundary
		B += sampler.E[index]
  end
  return abs(B)
end

# compute std of estimator using the indices in the index set
# TODO if we store the current index set in the sampler, we do not need to pass the index set
function compute_stochastic_error{d,N}(sampler, failProb, indexset::Set{Index{d,Vector{N}}})
  V = zeros(Float64,sampler.Z)
	for index::Index{d,Vector{N}} in indexset
		V += sampler.Vest[index]
  end
  return sqrt(2)*erfcinv(failProb)*sqrt(sum(V))
end

# compute mean of the estimator
# TODO same comment as above
function mean{d,N}(sampler, indexset::Set{Index{d,Vector{N}}})
  mu = zeros(Float64,sampler.Z)
	for index::Index{d,Vector{N}} in indexset
		mu += abs(sampler.E[index]) # TODO here we need a better estimate... abs or not? abs is consistent with regression model
  end
	return maximum(mu)
end

# compute std of the estimator
# TODO same comment as above
function std{d,N}(sampler, indexset::Set{Index{d,Vector{N}}})
  sigma = zeros(Float64,sampler.Z)
	for index::Index{d,Vector{N}} in indexset
		sigma += sampler.V[index]
  end
	return maximum(sigma)
end

# do regression on the indices provided based on the information already present in neighbours
# if the index falls inside the 3 x 3 x ... cube, do multilineair interpolation
# else, do repeated lineair interpolation and take the average
# linear regression of E / Vf / W for all given indices
function do_regression{d,N}(sampler, indices::Set{Index{d,Vector{N}}}) 
	Ds = [sampler.E, sampler.Vf, sampler.W]
	for ds in 1:length(Ds)
		corr = ds == 3 ? 1 : -1
		@debug println("----------------------------------------")				
		@debug which = ds == 1 ? "E" : ds == 2 ? "Vf" : "W"
		@debug println("performing regression on $(which)")
		@debug println("----------------------------------------")
		D = Ds[ds]
		# multi-linear regression of E / Vf / W for all given indices	
		Dcopy = copy(D) # make a copy of the dict...
		delete!(Dcopy, Index(zeros(N,d))) # ...and remove element 0
		@debug print("using ")
		@debug prettyprint(Set(keys(Dcopy)))
		(X,xi) = leastSquaresFit(Dcopy, 0) # multi-linear regression
		@debug println("results of multi-linear regression:")
		@debug println("  xi = $(xi)")
		@debug println("  X = $(X)")
		for index in indices # now use the regression model to get the value at the unknown indices
			@debug print_with_color(:cyan,"*** regression on index $(index.indices) ***\n")
			# or do repeated lineair regression if possible
			m = Float64[]
			for i in 1:d
							@debug print("checking if dimension $(i) is suitable for 1d regression... ")
				if index[i] > 2 # this index is suitable for repeated 1d regression
					@debug println("yes")
					# make a dict that contains the indices
					dict = typeof(D)()
					for j in 1:index[i]-1
						idx = copy(index)
						idx[i] = j
						@assert haskey(D,idx)
						dict[idx] = D[idx]
					end
					@debug print("using ")
					@debug prettyprint(Set(keys(dict)))
					# regress and store the result in m
					(Xacc,xiacc) = leastSquaresFit(dict,i)
					@debug println("results of 1d-linear regression:")
					@debug println("  xi = $(xiacc)")
					@debug println("  X = $(Xacc)")
					push!(m,Xacc*prod(2.^(corr*xiacc.*index)))
					@debug println("I guess the value from dimension $(i) is $(m[end])")
				else
					@debug println("no")
				end # if
			end # for
			val = isempty(m) ? X*prod(2.^(corr*xi.*index)) : mean(m)
			@debug println("regressed the value $(val) on index $(index.indices)")
			D[index] = eltype(values(D)) <: AbstractVector ? [val] : val
		end # for
	end # for
end # function

function t()

  # correction when doing continuation
  if continuate
    (A,α,B,β,Bacc,βacc,C,γ) = estimateProblemParameters(sampler)  # ! fits only max E_\ell, V_\ell
    Vftilde = bayesianUpdateVariance(sampler,A,α,B,β,E,Vf)
    Vtilde = bayesianUpdateVariance(sampler,A,α,Bacc,βacc,E,V)
    for index::Index{d,Vector{N}} in sort(indices) # use model fit on all indices
#						println("<<< WARNING >>> I ONLY USE A FIT ON THE BOUNDARY INDICES <<< END WARNING >>>")
			if !haskey(sampler.samples,index)
			E[index] = [A*prod(2.^(-α.*index))]
      Vf[index] = haskey(Vftilde,index) ? [Vftilde[index]] : [B*prod(2.^(-β.*index))]
      V[index] = haskey(Vtilde,index) ? [Vtilde[index]] : [Bacc*prod(2.^(-βacc.*index))]
      W[index] = C*prod(2.^(γ.*index))
      Vest[index] = haskey(sampler.samples,index) ? V[index]/size(sampler.samples[index],dir) : [zero(T)]
      if typeof(sampler.indexSet) == AD && haskey(profit,index)
	      profit[index] = maximum(abs(E[index])./sqrt(Vf[index].*W[index]))
      end
			end
    end
  end
	println("variance at index 0 is $(Vf[Index(zeros(N,d))])")
	# correct level 0
	#ndex = Index(zeros(N,d))
	#E[ndex] = squeeze(mean(sampler.samples[ndex],(1,2)),(1,2))
  #nsamples = size(sampler.samples[ndex],dir)
  #W[ndex] =  sampler.useTime ? sampler.times[ndex]/nsamples : prod(2.^(sampler.γ.*ndex))
  #Vest[ndex] = V[ndex]/nsamples
  #if typeof(sampler.indexSet) <: AD && haskey(profit,ndex) # compute gains when adaptive
	#  profit[ndex] = maximum(abs(E[ndex])./sqrt(Vf[ndex].*W[ndex]))
  #end
end

# estimate problem parameters A, α, B, β, Bacc, βacc, C, γ and Vtilde in the continuation problem
function estimateProblemParameters{d}(sampler::Sampler{d})

  # for now, let A, α, B, β, Bacc, βacc, C and γ follow from a fit through the available E_\ell, V_\ell and W_\ell
  E = Dict{Index{d,Vector{Int64}},Float64}()
  Vf = Dict{Index{d,Vector{Int64}},Float64}()
  V = Dict{Index{d,Vector{Int64}},Float64}()
  W = Dict{Index{d,Vector{Int64}},Float64}()
  for index in keys(sampler.samples)
    E[index] = maximum(squeeze(mean(sampler.samples[index],(1,2)),(1,2)))
    Vf[index] = maximum(squeeze(mean(var(sampler.samples[index],2),1),(1,2)))
    V[index] = maximum(squeeze(var(mean(sampler.samples[index],1),2),(1,2)))
    W[index] = sampler.useTime ? sampler.times[index] : prod(2.^(index.*sampler.γ))
  end

  # least-squares fit
  (A,α) = leastSquaresFit(E)
  (B,β) = leastSquaresFit(Vf)
  (Bacc,βacc) = leastSquaresFit(V)
  (C,γ) = leastSquaresFit(W)

  return (A::Float64, α::Vector{Float64}, B::Float64, β::Vector{Float64}, Bacc::Float64, βacc::Vector{Float64}, C::Float64, γ::Vector{Float64})
end

# Bayesian update of variances
function bayesianUpdateVariance{d,N1<:Integer,T<:AbstractFloat}(sampler::Sampler,A::T,α::Vector{T},B::T,β::Vector{T}, E::Dict{Index{d,Vector{N1}},Vector{T}}, V::Dict{Index{d,Vector{N1}},Vector{T}})
  dir = isa(sampler.numberGenerator,MCgenerator) ? 2 : 1
  Vtilde = Dict{Index{d,Vector{N1}},T}()

	for index in intersect(keys(E),keys(sampler.samples))
    N = size(sampler.samples[index],dir)
    μ = A*prod(2.^(-α.*index))
    λ = 1/B*prod(2.^(β.*index))
    Γ3 = sampler.k[2]*λ + N/2
    Γ4 = sampler.k[2] + 0.5*(N-1)*maximum(V[index]) + sampler.k[1]*N*(maximum(E[index])-μ)^2/(2*(sampler.k[1]+N))
    Vtilde[index] = Γ4/Γ3
  end
  zero_idx = Index(zeros(N1,d))
  Vtilde[zero_idx] = maximum(V[zero_idx]) # ! correction for level / index 0

  return Vtilde::Dict{Index{d,Vector{N1}},T}
end

# do a least-squares fit of the indices and data provided in the dict
# the optional argument d specifies the number of dimensions used in the interpolation
# if dim = 0, we do a d-lineair interpolation, else we do 1-d interpolation in the dim-th dimension 
function leastSquaresFit{d,N<:Integer,T}(dict::Dict{Index{d,Vector{N}},T}, dim::N)
	@assert dim >= 0 && dim <= d
  keyz = sort(Set(collect(keys(dict))))
	effdim = dim == 0 ? d : 1
	r = dim == 0 ? (1:d) : dim
  m = length(keyz)
  data = zeros(m,effdim+1)
  rhs = zeros(m)
  for i in 1:m # ! correction for level / index 0
		data[i,1:effdim] = keyz[i].indices[r]
    data[i,effdim+1] = 1
		rhs[i] = log2(abs(maximum(dict[keyz[i]][1])))
  end
	coeffs = data\rhs # least-squares system
  X = 2^coeffs[effdim+1]
  ξ = abs(coeffs[1:effdim])
	
	# X is constant, xi are d-dimensional rates			
  return (X,ξ)
end


# sample nbOfSamples from the Sampler at the given index
function sample{N<:Integer,I<:Index}(sampler::Sampler, nbOfSamples::N, index::I)

  # print some info to screen
  itype = ndims(sampler.indexSet) == 1 ? "level" : "index"
  idcs = ndims(sampler) == 1 ? 1 : 1:ndims(sampler)
  if isa(sampler.numberGenerator, MCgenerator)
      !sampler.showInfo || print("taking $(nbOfSamples) samples at "*itype*" $(index.indices[idcs])... ")
  else
      !sampler.showInfo || print("taking $(nshifts(sampler.numberGenerator))x$(nbOfSamples) samples at "*itype*" $(index.indices[idcs])... ")
  end

  # ask state of generator
  state = sampler.generatorStates
  if !haskey(state,index)
    state[index] = 0
  end
  
	# total number of workers to be used
  p = max(1,get(sampler.procMap,index,Void) - 1) # nworkers() = nprocs() - 1
	wp = WorkerPool(collect(1:p)+1)	# workerpool

	# parallel execution
	t = @elapsed r = pmap(wp, (i)->take_a_sample(index,state[index]+i,sampler), 1:nbOfSamples)

	# fetch results
  comb_samples = reduce(vcat,r)::Array{Float64,4}
  @assert size(comb_samples,1) == nbOfSamples
  samples = comb_samples[:,:,:,1]
  samples0 = sampler.storeSamples0 ? comb_samples[:,:,:,2] : zeros(0,0,0)

  # set state of generator
  state[index] += nbOfSamples

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
	println("done")

	# update dicts
	print("updating dicts... ")
  dir = isa(sampler.numberGenerator,MCgenerator) ? 2 : 1
	sampler.E[index] = squeeze(mean(sampler.samples[index],(1,2)),(1,2))
  sampler.Vf[index] = squeeze(mean(var(sampler.samples[index],2),1),(1,2))
  sampler.V[index] = squeeze(var(mean(sampler.samples[index],1),2),(1,2))
	n = size(sampler.samples[index],dir) # number of samples taken
  sampler.Vest[index] = sampler.V[index]/n
  sampler.Wst[index] =  prod(2.^(sampler.γ.*index))
	sampler.T[index] = t + get(sampler.T,index,0.0) # cumulative time 
  sampler.W[index] =  sampler.useTime ? sampler.T[index]/n : sampler.Wst[index]
	sampler.P[index] = maximum(abs(sampler.E[index])./sqrt(sampler.Vf[index].*sampler.W[index]))

	println("done")
  return Void
end

# take a single sample from the sampler, function to be run in parallel using pmap
function take_a_sample{N<:Integer,I<:Index}(index::I, i::N, sampler::Sampler)

  # generate "random" numbers
  XI = getPoint(sampler.numberGenerator,i)

  # compute difference
	mySample = zeros(Float64,1,nshifts(sampler.numberGenerator),sampler.Z,2)::Array{Float64,4}
	for q = 1:nshifts(sampler.numberGenerator)
    xi = XI[:,q]::Vector{Float64}
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
	
end
