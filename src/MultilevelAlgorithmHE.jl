# simulate(sampler, absTOL)
simulate{T<:AbstractFloat}(sampler::Sampler, absTOL::T) = simulate(sampler, absTOL=absTOL)

# simulate(sampler, absTOL, failProb)
simulate{T<:AbstractFloat}(sampler::Sampler, absTOL::T, failProb::T) = simulate(sampler, absTOL=absTOL, failProb=failProb)

# simulate(sampler, absTOL=..., relTOL=..., failProb=...)
function simulate{T<:AbstractFloat}(sampler::Sampler; absTOL::T=Inf, relTOL::T=Inf, failProb::T=0.1)

  # checks on inputs
  if absTOL < 0 || relTOL < 0
    error("supplied tolerances (absolute and/or relative) must be positive!")
  elseif failProb < 0 || failProb > 1
    error("failure probability failProb must be between 0 and 1!")
  elseif isinf(absTOL) && isinf(relTOL)
    error("supply a relative or absolute tolerance!")
  end

  t = Vector{T}()

  N = Int64
  old = Set{Index{3,Vector{N}}}()
  profit = Dict{Index{3,Vector{N}},T}()
  profit[Index(zeros(N,3))] = 0.
  L = 0

  # run mimc or cmimc
  if sampler.continuate
    absTOL_ = 2^(sampler.nTOL-1)*absTOL
    relTOL_ = 2^(sampler.nTOL-1)*relTOL
		relTOL_ = 0.16
    sampler.continuate = false # run initial hierarchy
    tic()
    #old, profit, L = mimc(sampler,absTOL_,relTOL_,failProb,old,profit,L)
    minL = mimc(sampler,absTOL_,relTOL_,failProb,old,profit,L)
    push!(t, toc())
    sampler.continuate = true
    
		#for i in 1:sampler.nTOL-1
		#tols = [0.0800    0.0400    0.0200 0.016  0.012  0.0100  0.009 0.008 0.0050]
		tols=[0.0800 0.0400  0.03 0.024  0.0200 0.016  0.012  0.0100  0.009 0.008 0.0060]
		#tols = [0.012,0.01,0.008,0.006,0.004,0.002]
		for i in 1:length(tols)
			#absTOL_ = 2^(sampler.nTOL-1-i)*absTOL
      #relTOL_ = 2^(sampler.nTOL-1-i)*relTOL
      tic()
			relTOL_ = tols[i]
absTOL_ = Inf
      #old, profit, L = mimc(sampler,absTOL_,relTOL_,failProb,old,profit,L)
      minL = mimc(sampler,absTOL_,relTOL_,failProb,old,profit,minL)
      
			endtime = toc()
			println("ELAPSED IS $(endtime)")
			
			if endtime > 24*3600
							println("TIME EXCEEDED MAX TIME -- QUITTING...")
							println("TIME EXCEEDED MAX TIME -- QUITTING...")
							println("TIME EXCEEDED MAX TIME -- QUITTING...")
							println("TIME EXCEEDED MAX TIME -- QUITTING...")
							println("TIME EXCEEDED MAX TIME -- QUITTING...")
							println("TIME EXCEEDED MAX TIME -- QUITTING...")
				break
			end
			
			push!(t, endtime)
      #save(sampler)
    end
  else
      push!(t, @elapsed mimc(sampler,absTOL,relTOL,failProb,old,profit,L)) # plain mimc
  end

  return cumsum(t) # return timing results

end

# mimc simulation
function mimc{d,T<:AbstractFloat,N<:Integer}(sampler::Sampler{d}, absTOL::T, relTOL::T, failProb::T, old::Set{Index{d,Vector{N}}}, profit::Dict{Index{d,Vector{N}},T}, minL::N)

  # print some info to screen
  if ( sampler.showInfo )
    @printf(sampler.ioStream,"%s","--------------------------------------------------------------------------------\n")
    @printf(sampler.ioStream,"%s","*** MultilevelEstimators.jl @$(now())\n")
    @printf(sampler.ioStream,"%s","*** Simulating $(sampler.sampleFunction)\n")
    isCont = sampler.continuate ? "" : "not"
    idxSet = ( isa(sampler.indexSet,ML) || isa(sampler.indexSet,SL) || isa(sampler.indexSet,AD) ) ? "$(sampler.indexSet)" : "$(sampler.indexSet)"[1:34]
    @printf(sampler.ioStream,"%s","*** Using a $(idxSet), $isCont continuating \n")
    @printf(sampler.ioStream,"%s",
    @sprintf("*** absTOL = %0.3e / relTOL = %0.3e (failure probability of %0.2f)\n",absTOL,relTOL,failProb))
    @printf(sampler.ioStream,"%s","--------------------------------------------------------------------------------\n")
  end

  # loop variables
  converged = false
  #
  L  = 0#typeof(sampler.indexSet) <: AD ? length(sampler.samples) : 0
  #

  # aliases
#  N = Int64
  λ = sampler.numberGenerator.λ
  q = nshifts(sampler.numberGenerator)
  dir = isa(sampler.numberGenerator,MCgenerator) ? 2 : 1

  # variable definition
  F = zeros(T,sampler.Z); G = zeros(T,sampler.Z)
  E = Dict{Index{d,Vector{N}},Vector{T}}() # dict for expected values
  V = Dict{Index{d,Vector{N}},Vector{T}}() # dict for variances
  Vf = Dict{Index{d,Vector{N}},Vector{T}}() # dict for variance of the function
  W = Dict{Index{d,Vector{N}},T}() # dict for costs
  S = Dict{Index{d,Vector{N}},N}() # dict for optimal number of samples
  Vest = Dict{Index{d,Vector{N}},Vector{T}}() # dict for variance of estimator
  C = sqrt(2)*erfcinv(failProb)
  oldindexset = Set{Index{d,Vector{N}}}()
  #
  old = Set{Index{d,Vector{N}}}()
  profit = Dict{Index{d,Vector{N}},T}()
  profit[Index(zeros(N,d))] = 0.
  #
  dr = zeros(0,d)


  while !converged
    # empty dicts
    E = Dict{Index{d,Vector{N}},Vector{T}}() # dict for expected values
    V = Dict{Index{d,Vector{N}},Vector{T}}() # dict for variances
    Vf = Dict{Index{d,Vector{N}},Vector{T}}() # dict for variance of the function
    W = Dict{Index{d,Vector{N}},T}() # dict for costs
    S = Dict{Index{d,Vector{N}},N}() # dict for optimal number of samples
    Vest = Dict{Index{d,Vector{N}},Vector{T}}() # dict for variance of estimator




    # get index set in current iteration
    if typeof(sampler.indexSet) <: AD
	#for index in old
	#	profit[index] = 0.
	#end
      if L<=3 #&& isempty(keys(sampler.samples))
	newindexset = getIndexSet(TD(d),L)
	newindexsettorem = getIndexSet(TD(d),L-1)
	for index in newindexset
		profit[index]=0.
	  dr = vcat(dr,(index.indices)')
	  writedlm("data/indexsetL$(L).txt",dr)
	end
	for index in newindexsettorem

		delete!(profit,index)
		push!(old,index)
	end
      else
        # find index with largest gain (sparse addition)
	i = collect(keys(profit))[indmax(collect(values(profit)))]
	println("*****************************************************")
	println("*****************************************************")
	
	println("                   ::::PROFITS::::")	
	println(profit)	
	println("*****************************************************")
	println("*****************************************************")
	println("===== adding $(i) =====")
	delete!(profit,i) # delete from profit set
        push!(old,i) # add to old set and check admissables
        for p in 1:d 
          j = copy(i)::Index{d,Vector{N}}
          j[p] += 1
					if isAdmissable(union(old, Set(keys(profit))),j)
		println("!!! MIND !!! hacked shape of index set: max is (10,9,9)")
		if !(j[1] > 10 || j[2] > 9 || j[3] > 9)
		  profit[j] = 0.
	  dr = vcat(dr,(j.indices)')
	  writedlm("data/indexsetL$(L).txt",dr)
	  	end
          end
        end
        newindexset = union(old, Set(keys(profit)))::Set{Index{d,Vector{N}}}
      end
    else
      newindexset = getIndexSet(sampler.indexSet,L)::Set{Index{d,Vector{N}}}
    end

#
# ITERATE OVER ALL INDICES TO DETERMINE OLD
# i.e. move additional indices from profit to old
if typeof(sampler.indexSet)<:AD
				active = Set(keys(profit))
				for index in active
								for j in 1:d
												e = zeros(N,d)
												e[j] = 1
												neig = Index(index.indices+e)
												if in(neig,newindexset)
																delete!(profit,index)
																push!(old,index)
																println("*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/")
																println("$(index) removed!")
																println("*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/")
												end
								end
				end
end
	
    # get boundary in current iteration
    boundary = Set{Index{d,Vector{N}}}()::Set{Index{d,Vector{N}}} # empty set
    if typeof(sampler.indexSet) <: ML
      push!(boundary,Index(L*ones(N,d)))
    elseif typeof(sampler.indexSet) <: AD
      #boundary = union(boundary,getBoundary(newindexset)::Set{Index{d,Vector{N}}})
      boundary = Set(keys(profit))
 #     println("============================================================")
 #     println("============================================================")
 #     println("            BOUNDARY")
 #     println(boundary)
 #     println("============================================================")
 #     println("============================================================")
      
    else
      boundary = union(boundary,getBoundary(getIndexSet(sampler.indexSet, L))::Set{Index{d,Vector{N}}})
    end
  #     println("============================================================")
  #    println("============================================================")
  #    println("            BOUNDARY")
  #    println(boundary)
  #    println("============================================================")
  #    println("============================================================")
    
    # new indices to add
    indicesToAdd = setdiff( newindexset, oldindexset )


   # println("indicesToAdd $(indicesToAdd)")

    # print some info to the screen
    !sampler.showInfo || print_with_color(:blue, sampler.ioStream, "*** currently running at L = $(L)...\n")

    # take initial number of samples
    if !sampler.continuate
      for index::Index{d,Vector{N}} in indicesToAdd
        if !haskey(sampler.times,index)
          time = @elapsed sample(sampler, sampler.Nstar, index)
          setindex!(sampler.times, time, index) # cumulative time
#	  dr = vcat(dr,(index.indices)')
#	  writedlm("data/indexsetL$(L).txt",dr')
  	end
      end
    end
    updateDicts!(sampler,newindexset,false,E,Vf,V,W,Vest,profit)

    # determine actual absolute tolerance
    #if L=0
	    
    minTOL = min(absTOL,relTOL*maximum(sum(values(E))))
    #minTOL = min(absTOL,relTOL*maximum(sum(values(E))))
#else

  #  minTOL = min(absTOL,relTOL*mye)
 #   end
    

 ######
 ######
 ######
 ###### DO NOT CHANGE THIS LINE OR ERROR 
 ######
 ######
    updateDicts!(sampler,newindexset,sampler.continuate,E,Vf,V,W,Vest,profit)
#    updateDicts!(sampler,newindexset,false,E,Vf,V,W,Vest,profit)
    # estimate initial bias and splitting
    Gf = zeros(T,sampler.Z)
    for index::Index{d,Vector{N}} in boundary
						Gf += abs(E[index]) #abs(E[index])
    end
    G = maximum(abs(Gf))
    !sampler.showInfo || print_with_color(:green, sampler.ioStream, @sprintf("   *** initial bias = %0.6e \n", G ) )
    splitting = ( G < minTOL/2 ) ? 1 - G/minTOL : 0.5
    !sampler.showInfo || print_with_color(:green, sampler.ioStream, @sprintf("   *** splitting = %0.6e \n", splitting ) )

    updateDicts!(sampler,newindexset,sampler.continuate,E,Vf,V,W,Vest,profit)






  #if ( sampler.showInfo )
  #  @printf(sampler.ioStream,"%s","--------------------------------------------------------------------------------\n")
  #  itype = ndims(sampler.indexSet) == 1 ? "level" : "index"
  #  @printf(sampler.ioStream,"%s","  "*itype*"       E              V               N               W              \n")
  #  @printf(sampler.ioStream,"%s","--------------------------------------------------------------------------------\n")
  #  for index::Index{d,Vector{N}} in sort(newindexset)
  #    str = ("  $(index.indices)            "[1:13])
  #    str *= @sprintf("%12.5e",maximum(E[index]))
  #    str *= @sprintf("    %0.6e",maximum(Vf[index]))
  #    str *= sampler.continuate ? "" : @sprintf("    %d               ",prod(size(sampler.samples[index])[1:2]))[1:16]
  #    str *= @sprintf("    %0.6e",W[index])
  #    str *= @sprintf("    %0.6e\n",maximum(abs(E[index])./sqrt(Vf[index].*W[index])))
  #    @printf(sampler.ioStream,"%s",str)
  #  end
  #end










    # calculate optimal number of samples
    mySum = zeros(T,sampler.Z)
    for index::Index{d,Vector{N}} in newindexset
      mySum += ( ( Vf[index].*(W[index]*ones(sampler.Z)).^(2*λ) ).^(1/(2*λ+1) ) )
    end
    for index::Index{d,Vector{N}} in newindexset
      Nopt = ( ( C/(splitting*minTOL) )^2 * 1/q * (Vf[index]./W[index]).^(2*λ/(2*λ+1)) .* mySum ).^(1/(2*λ))
      S[index] = ceil(N,max(3.,maximum(Nopt)))
    end
   # println("-----------------")
	#println(S)
   # println("-----------------")

    # take additional samples at each level
    for index::Index{d,Vector{N}} in newindexset
      samplesToTake = haskey(sampler.samples,index) ? max( 0, S[index]-size(sampler.samples[index],dir) ) : S[index]
      if samplesToTake > 0
        time = @elapsed sample(sampler, samplesToTake, index )
        sampler.times[index] = haskey(sampler.times,index) ? sampler.times[index] + time : time # cumulative time
      end
    end
    #updateDicts!(sampler,newindexset,sampler.continuate,E,Vf,V,W,Vest,profit)
####
###
###
###
#### THIS LINE
###
###
###
####
    updateDicts!(sampler,newindexset,false,E,Vf,V,W,Vest,profit)
    # safety
    Gf = zeros(T,sampler.Z)
    for index::Index{d,Vector{N}} in boundary
						Gf += abs(E[index]) #abs(E[index])
    end
    G = maximum(abs(Gf))
    updateDicts!(sampler,newindexset,sampler.continuate,E,Vf,V,W,Vest,profit)
    F = C*sqrt(sum(values(Vest)))
    #G = zero(T)
    #for index::Index{d,Vector{N}} in boundary
    #  G += maximum(abs(E[index]))
    #end
    splitting = ( G < minTOL/2 ) ? 1 - G/minTOL : 0.5
    if sampler.safety
      while maximum(F) > splitting*minTOL
        # double number of samples on index with max ratio Vest/W
        maxratio = zero(T)
        maxindex = Index(zeros(N,d))::Index{d,Vector{N}}
        for index::Index{d,Vector{N}} in newindexset
          ratio = maximum(Vest[index])/W[index]
          if ratio > maxratio
            maxratio = ratio
            maxindex = index::Index{d,Vector{N}}
          end
        end
        n = size(sampler.samples[maxindex],dir)
        time = @elapsed sample(sampler, nextpow2(n+1)-n, maxindex ) # round to nearest power of two
        sampler.times[maxindex] += time # cumulative time
        updateDicts!(sampler,newindexset,sampler.continuate,E,Vf,V,W,Vest,profit)
##########################################################
# referee 1 addition
    minTOL = min(absTOL,relTOL*maximum(sum(values(E))))
# calculate optimal number of samples
    mySum = zeros(T,sampler.Z)
    for index::Index{d,Vector{N}} in newindexset
      mySum += ( ( Vf[index].*(W[index]*ones(sampler.Z)).^(2*λ) ).^(1/(2*λ+1) ) )
    end
    for index::Index{d,Vector{N}} in newindexset
      Nopt = ( ( C/(splitting*minTOL) )^2 * 1/q * (Vf[index]./W[index]).^(2*λ/(2*λ+1)) .* mySum ).^(1/(2*λ))
      S[index] = ceil(N,max(3.,maximum(Nopt)))
    end

    # take additional samples at each level
    for index::Index{d,Vector{N}} in newindexset
      samplesToTake = haskey(sampler.samples,index) ? max( 0, S[index]-size(sampler.samples[index],dir) ) : S[index]
      if samplesToTake > 0
        time = @elapsed sample(sampler, samplesToTake, index )
        sampler.times[index] = haskey(sampler.times,index) ? sampler.times[index] + time : time # cumulative time
      end
    end
    updateDicts!(sampler,newindexset,sampler.continuate,E,Vf,V,W,Vest,profit)
##########################################################
	
	# recompute statistical error and splitting
        F = C*sqrt(sum(values(Vest)))
     Gf = zeros(T,sampler.Z)
    for index::Index{d,Vector{N}} in boundary
						Gf += abs(E[index]) #abs(E[index])
    end
    G = maximum(abs(Gf))
 #G = zero(T)
  #      for index::Index{d,Vector{N}} in boundary
   #       G += maximum(abs(E[index]))
    #    end
        splitting = ( G < minTOL/2 ) ? 1 - G/minTOL : 0.5
      end
    end

    # check for convergence
    updateDicts!(sampler,Set(collect(keys(sampler.samples))),false,E,Vf,V,W,Vest,profit) # store true values in Dicts
		if ( L > 1 ) && ( L >= minL )
      converged = ( maximum(F + G)::T < minTOL )
      println("minTOL is $minTOL")
      println("we have error $(maximum(F+G))")
      # print some info to the screen
      !sampler.showInfo || print_with_color(:magenta, sampler.ioStream, 
        @sprintf("*** error estimate is %0.6e + %0.6e = %0.6e \n", maximum(F), maximum(G), maximum(F+G) ) )
      μ = maximum(sum(values(E)))
      σ = sqrt(maximum(sum(values(V))))
      !sampler.showInfo || ( !converged || print_with_color(:magenta, sampler.ioStream, 
        @sprintf("*** result is %0.6e ±(%0.6e, %0.6e, %0.6e) \n", μ, σ, 2*σ, 3*σ ) ) )
    end

#    mye = maximum(sum(values(E))


    # print warning if no convergence
    if !converged
      if L == sampler.maxL
        warn("maximum level reached and no convergence yet, sorry! :(\n")
        converged = true
      end
    end

    # update L
    L += 1
    oldindexset = newindexset
  ######3end

  #
  # PRINT OVERVIEW TABLE
  #
  if ( sampler.showInfo )
    @printf(sampler.ioStream,"%s","--------------------------------------------------------------------------------\n")
    itype = ndims(sampler.indexSet) == 1 ? "level" : "index"
    @printf(sampler.ioStream,"%s","  "*itype*"       E              V               N               W              \n")
    @printf(sampler.ioStream,"%s","--------------------------------------------------------------------------------\n")
    for index::Index{d,Vector{N}} in sort(newindexset)
      str = ("  $(index.indices)            "[1:13])
      str *= @sprintf("%12.5e",maximum(E[index]))
      str *= @sprintf("    %0.6e",maximum(Vf[index]))
      str *= @sprintf("    %d               ",prod(size(sampler.samples[index])[1:2]))[1:16]
      str *= @sprintf("    %0.6e",W[index])
      str *= @sprintf("    %0.6e\n",maximum(abs(E[index])./sqrt(Vf[index].*W[index])))
      @printf(sampler.ioStream,"%s",str)
    end
  end


  end

  return L-1
end

# refactored function that updates E/Vf/V/W/Vest/profit from current sampler
function updateDicts!{d,N,T}(sampler::Sampler, indices::Set{Index{d,Vector{N}}}, continuate::Bool,
  E::Dict{Index{d,Vector{N}},Vector{T}}, Vf::Dict{Index{d,Vector{N}},Vector{T}}, 
    V::Dict{Index{d,Vector{N}},Vector{T}}, W::Dict{Index{d,Vector{N}},T}, 
      Vest::Dict{Index{d,Vector{N}},Vector{T}}, profit::Dict{Index{d,Vector{N}},T})

  dir = isa(sampler.numberGenerator,MCgenerator) ? 2 : 1

  # update all
  for index::Index{d,Vector{N}} in sort(indices)
    if haskey(sampler.samples,index)
      E[index] = squeeze(mean(sampler.samples[index],(1,2)),(1,2))
      Vf[index] = squeeze(mean(var(sampler.samples[index],2),1),(1,2))
      V[index] = squeeze(var(mean(sampler.samples[index],1),2),(1,2))
      nsamples = size(sampler.samples[index],dir)
      W[index] =  sampler.useTime ? sampler.times[index]/nsamples : prod(2.^(sampler.γ.*index))
      Vest[index] = V[index]/nsamples
      if typeof(sampler.indexSet) <: AD && haskey(profit,index) # compute gains when adaptive
	      #profit[index] = min(1,maximum(abs(E[index]./E[Index(0,0,0)])./sqrt(Vf[index]./Vf[Index(0,0,0)].*W[index]./W[Index(0,0,0)])))
	      profit[index] = maximum(abs(E[index])./sqrt(Vf[index].*W[index]))
      end
    end
  end

  # correction when doing continuation
  if continuate
    (A,α,B,β,Bacc,βacc,C,γ) = estimateProblemParameters(sampler)  # ! fits only max E_\ell, V_\ell
    Vftilde = bayesianUpdateVariance(sampler,A,α,B,β,E,Vf)
    Vtilde = bayesianUpdateVariance(sampler,A,α,Bacc,βacc,E,V)
    for index::Index{d,Vector{N}} in sort(indices) # use model fit on all indices
      E[index] = [A*prod(2.^(-α.*index))]
      Vf[index] = haskey(Vftilde,index) ? [Vftilde[index]] : [B*prod(2.^(-β.*index))]
      V[index] = haskey(Vtilde,index) ? [Vtilde[index]] : [Bacc*prod(2.^(-βacc.*index))]
      W[index] = C*prod(2.^(γ.*index))
      Vest[index] = haskey(sampler.samples,index) ? V[index]/size(sampler.samples[index],dir) : [zero(T)]
      if typeof(sampler.indexSet) == AD && haskey(profit,index)
	      profit[index] = maximum(abs(E[index])./sqrt(Vf[index].*W[index]))
#	      profit[index] = min(1,maximum(abs(E[index]./E[Index(0,0,0)])./sqrt(Vf[index]./Vf[Index(0,0,0)].*W[index]./W[Index(0,0,0)])))
      end
    end
  end




#correct level 0
  ndex = Index([0,0,0])
E[ndex] = squeeze(mean(sampler.samples[ndex],(1,2)),(1,2))
      #Vf[index] = squeeze(mean(var(sampler.samples[index],2),1),(1,2))
      #V[index] = squeeze(var(mean(sampler.samples[index],1),2),(1,2))
      nsamples = size(sampler.samples[ndex],dir)
      W[ndex] =  sampler.useTime ? sampler.times[ndex]/nsamples : prod(2.^(sampler.γ.*ndex))
      Vest[ndex] = V[ndex]/nsamples
      if typeof(sampler.indexSet) <: AD && haskey(profit,ndex) # compute gains when adaptive
	      #profit[index] = min(1,maximum(abs(E[index]./E[Index(0,0,0)])./sqrt(Vf[index]./Vf[Index(0,0,0)].*W[index]./W[Index(0,0,0)])))
	      profit[ndex] = maximum(abs(E[ndex])./sqrt(Vf[ndex].*W[ndex]))
      end
    











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
  println("alpha = $(α)")
  (B,β) = leastSquaresFit(Vf)
  println("beta = $(β)")
  (Bacc,βacc) = leastSquaresFit(V)
  println("beta_acc = $(βacc)")
  (C,γ) = leastSquaresFit(W)
  println("gamma = $(γ)")

  # println("I have fitted A = $(A) and α= $(α)")
  # println("--------------------------------------------------------------------------------")
  # itype = ndims(sampler.indexSet) == 1 ? "level" : "index"
  # println("  "*itype*"       E              E_model         error                  ")
  # println("--------------------------------------------------------------------------------")
  # for index in sort(Set(collect(keys(E))))
  #   str = ("  $(index.indices)            "[1:13])
  #   str *= @sprintf("%12.5e",maximum(abs(E[index])))
  #   str *= @sprintf("    %0.6e",A*prod(2.^(-α.*index.indices)))
  #   str *= @sprintf("    %0.6e",abs(A*prod(2.^(-α.*index.indices))-maximum(abs(E[index])))/(maximum(abs(E[index]))))
  #   println(str)
  # end

  # println("I have fitted B = $(B) and β= $(β)")
  # println("--------------------------------------------------------------------------------")
  # itype = ndims(sampler.indexSet) == 1 ? "level" : "index"
  # println("  "*itype*"       Vf             Vf_model        error                  ")
  # println("--------------------------------------------------------------------------------")
  # for index in sort(Set(collect(keys(Vf))))
  #   str = ("  $(index.indices)            "[1:13])
  #   str *= @sprintf("%12.5e",maximum(abs(Vf[index])))
  #   str *= @sprintf("    %0.6e",B*prod(2.^(-β.*index.indices)))
  #   str *= @sprintf("    %0.6e",abs(B*prod(2.^(-β.*index.indices))-maximum(abs(Vf[index])))/(maximum(abs(Vf[index]))))
  #   println(str)
  # end

  # println("I have fitted Bacc = $(Bacc) and βacc = $(βacc)")
  # println("--------------------------------------------------------------------------------")
  # itype = ndims(sampler.indexSet) == 1 ? "level" : "index"
  # println("  "*itype*"       V              V_model         error                  ")
  # println("--------------------------------------------------------------------------------")
  # for index in sort(Set(collect(keys(V))))
  #   str = ("  $(index.indices)            "[1:13])
  #   str *= @sprintf("%12.5e",maximum(abs(V[index])))
  #   str *= @sprintf("    %0.6e",Bacc*prod(2.^(-βacc.*index.indices)))
  #   str *= @sprintf("    %0.6e",abs(Bacc*prod(2.^(-βacc.*index.indices))-maximum(abs(V[index])))/(maximum(abs(V[index]))))
  #   println(str)
  # end

  # println("I have fitted C = $(C) and γ= $(γ)")
  # println("--------------------------------------------------------------------------------")
  # itype = ndims(sampler.indexSet) == 1 ? "level" : "index"
  # println("  "*itype*"       W              W_model         error                  ")
  # println("--------------------------------------------------------------------------------")
  # for index in sort(Set(collect(keys(W))))
  #   str = ("  $(index.indices)            "[1:13])
  #   str *= @sprintf("%12.5e",maximum(abs(W[index])))
  #   str *= @sprintf("    %0.6e",C*prod(2.^(γ.*index.indices)))
  #   str *= @sprintf("    %0.6e",abs(C*prod(2.^(γ.*index.indices))-maximum(abs(W[index])))/(maximum(abs(W[index]))))
  #   println(str)
  # end

  # # index set from which to derive the parameters
  # L = max(1,sampler.L-sampler.bayesianProperties.ell_0)
  # I = setdiff(getIndexSet(sampler.indexset,sampler.L),getIndexSet(sampler.indexset,L))
  # numerator = 0.
  # denominator = 0.
  # for i in I
  #   numerator += prod(2^(-(α-β).*i.indices))*N[i]*E[i]
  #   denominator += N[i]*prod(2^(-(2*α-β).*i.indices))
  # end
  # A_ = numerator/denominator

  # numerator = 0.
  # denominator = 0.
  # for i in I
  #   numerator += prod(2^(β.*i.indices))*      N[i]*E[i]
  #   denominator += N[i]*prod(2^(-(2*α-β).*i.indices))
  # end
  # B_ = numerator/denominator

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
function leastSquaresFit{d,N<:Integer,T<:AbstractFloat}(dict::Dict{Index{d,Vector{N}},T})
  keyz = sort(Set(collect(keys(dict))))
  m = length(keyz)
  data = zeros(m-1,d+1)
  rhs = zeros(m-1)
  for i in 1:m-1 # ! correction for level / index 0
    data[i,1:d] = keyz[i+1].indices
    data[i,d+1] = 1
    rhs[i] = log2(abs(dict[keyz[i+1]][1]))
  end
  coeffs = data\rhs # least-squares system
  X = 2^coeffs[d+1]
  ξ = abs(coeffs[1:d])

  return (X,ξ)
end
