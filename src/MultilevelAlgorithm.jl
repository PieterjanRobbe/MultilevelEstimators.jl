#
# MULTILEVELALGORITHM.JL
#

#	Collection of functions that implement the different algorithms, MLMC, MLQMC, MIMC, MIQMC...
# Main function is "simulate", see below.

# simulate(sampler, absTol)
simulate{T<:AbstractFloat}(sampler::Sampler, absTol::T) = simulate(sampler, [absTol])

# simulate(sampler, absTol, failProb)
simulate{T<:AbstractFloat}(sampler::Sampler, absTol::T, failProb::T) = simulate(sampler, [absTol], failProb=failProb)

# simulate(sampler, tol, failProb=..., is_relative=..., folder=...)
# simulate takes a vector tol as input that contains the tolerances that must be solved for
function simulate{T<:AbstractFloat,N}(sampler::Sampler, tol::Vector{T}, failProb::T=0.1, is_relative::Bool=false, folder::AbstractString=".", max_time::N=typemax(Int64))

  # checks on inputs
	if any(tol .< 0) || any(isinf(tol))
    error("supplied tolerances must be positive!")
  elseif failProb < 0 || failProb > 1
    error("failure probability failProb must be between 0 and 1!")
	elseif !isdir(folder)
		error("unknown folder $(folder)!")
  end # if

	# if user provided a set of tolerances, sort them
	sort!(tol,rev=true)

	# else, use "optimal" sequence from CMLMC paper
	if sampler.continuate && length(tol) == 1
		tol = 1.5.^(sampler.nTOL-(1:sampler.nTOL))*tol[1]
	end # if

	# preallocate vector that will contain wall clock and standard cost
	wctime = zeros(tol)
	stcost = zeros(tol)

  # run mimc 
	for i in 1:length(tol)
		# simulate
		delta_t = @elapsed mimc(sampler,tol[i],is_relative,failProb) 
		# compute wall clock time and standard cost
		wctime[i] =  delta_t
		for index in keys(sampler.samples)
			stcost[i] += length(sampler.samples[index])*sampler.Wst[index]
		end # for
		@debug begin
			if !isdir("data")
				mkdir("data")
			end
			if !isdir("data/indexsets")
				mkdir("data/indexsets")
			end
			print_with_color(:cyan,"ELAPSED IS $(delta_t)\n")
			print("writing wall clock times into $(folder)/data/wctime.txt... ")
			writedlm(folder*"/data/wctime.txt",cumsum(wctime))
			println("done")
			print("writing standard cost into $(folder)/data/stcost.txt... ")
			writedlm(folder*"/data/stcost.txt",cumsum(stcost))
			println("done")
		end # begin
		if cumsum(wctime)[end] > max_time
			print_with_color(:red,"wall clock time exceeded maximum time of $(max_time) secondes, aborting...\n")
			break
		end # if
	end # for

	return cumsum(wctime),cumsum(stcost)
end # function

# actual mimc simulation
function mimc{d,T<:AbstractFloat}(sampler::Sampler{d}, TOL::T, is_relative::Bool, failProb::T)

  # print some info to screen
  if ( sampler.showInfo )
    @printf(sampler.ioStream,"%s","--------------------------------------------------------------------------------\n")
    @printf(sampler.ioStream,"%s","*** MultilevelEstimators.jl @$(now())\n")
    @printf(sampler.ioStream,"%s","*** Simulating $(sampler.sampleFunction)\n")
    isCont = sampler.continuate ? "" : "not"
    idxSet = ( isa(sampler.indexSet,ML) || isa(sampler.indexSet,SL) || isa(sampler.indexSet,AD) ) ? "$(sampler.indexSet)" : "$(sampler.indexSet)"[1:34]
    @printf(sampler.ioStream,"%s","*** Using a $(idxSet), $isCont continuating \n")
    printTOL = is_relative ? "relTOL" : "absTOL"
    @printf(sampler.ioStream,"%s","*** $printTOL = "*@sprintf("%0.3e (failure probability of %0.2f)\n",TOL,failProb))
		@printf(sampler.ioStream,"%s","--------------------------------------------------------------------------------\n")
  end # if

  # loop variables
  converged = false
  L  = 0

	# when debugging, keep track of parameter L, bias, stochastic error, total error
	@debug Lcollection = zeros(T,sampler.maxL+1,3)

  # aliases
	# TODO might be incorporated into the sampler object
  N = Int64
  λ = sampler.numberGenerator.λ
  q = nshifts(sampler.numberGenerator)
  dir = isa(sampler.numberGenerator,MCgenerator) ? 2 : 1

  # variable definition
  oldindexset = Set{Index{d,Vector{N}}}()

	# variables used when adaptive
  old = Set{Index{d,Vector{N}}}() # old index set
  active = Dict{Index{d,Vector{N}},T}() # active set
  active[Index(zeros(N,d))] = 0.
  
  while !converged
    # print some info to the screen
    !sampler.showInfo || print_with_color(:blue, sampler.ioStream, "*** currently running at L = $(L)...\n")
    
		# empty dicts
    S = Dict{Index{d,Vector{N}},N}() # dict for optimal number of samples

    # get index set in current iteration
		if !(typeof(sampler.indexSet) <: AD)
			newindexset = getIndexSet(sampler.indexSet,L)::Set{Index{d,Vector{N}}}	
		elseif typeof(sampler.indexSet) <: AD && L <= 2
						newindexset = getIndexSet(TD(d),L)::Set{Index{d,Vector{N}}}
		else
			if L <= 2 # assume TD index set for inial estimates
				newindexset = getIndexSet(TD(d),L)
				newindexsetboundary = getBoundary(newindexset)
				for index in newindexsetboundary
					active[index]=0.
				end
				for index in setdiff(newindexset, newindexsetboundary)
					push!(old,index)
				end
      else
        # find index with largest gain
				i = collect(keys(active))[indmax(collect(values(profit)))]
				delete!(active,i) # delete from active set
        push!(old,i) # add to old set and check admissables
        for p in 1:d 
          j = copy(i)::Index{d,Vector{N}}
          j[p] += 1
          if isAdmissable(old,j) # add when admissible in old set
		  			active[j] = 0.
	  			end
        end
        newindexset = union(old, Set(keys(profit)))::Set{Index{d,Vector{N}}}
      end
    end
		# print index set
		@debug print("the new index set is an ")
		@debug prettyprint(newindexset)
		# store the index set in a file
		@debug begin
			# store all indices in a matrix M
			M = zeros(T, length(newindexset), d)
			cntr = 1
			for idx in sort(newindexset)
				M[cntr,1:d] = idx.indices
				cntr += 1
			end
			writedlm("data/indexsets/L$(L).txt",M)
		end

    # get boundary in current iteration
    boundary = Set{Index{d,Vector{N}}}()::Set{Index{d,Vector{N}}} # empty set
    if typeof(sampler.indexSet) <: ML
      push!(boundary,Index(L*ones(N,d)))
    elseif typeof(sampler.indexSet) <: AD
      boundary = Set(keys(profit))
    else
      boundary = union(boundary,getBoundary(getIndexSet(sampler.indexSet, L))::Set{Index{d,Vector{N}}})
    end
		@debug print("the boundary is an ")
		@debug prettyprint(boundary)
    
    # new indices to add
    indicesToAdd = setdiff( newindexset, oldindexset )
		@debug print("the new indices are an ")
		@debug prettyprint(indicesToAdd)

    # take initial number of samples at each new index when needed...
    if L < 3
			@debug println("now taking warm-up samples... ")
      for index::Index{d,Vector{N}} in indicesToAdd
        if !haskey(sampler.T,index) # when running repeatedly, might already have samples available
          sample(sampler, sampler.Nstar, index)
  			end
      end
			@debug println("done taking warm-up samples!")
		# ... or else do regression over the new indices to get estimates for variances
		else
			@debug println("now doing regression... ")
			do_regression(sampler,indicesToAdd)
			@debug println("regression ok!")
		end
		@debug inspect(sampler)

		# actual relative/absolute tolerance to solve for
    realTOL = is_relative ? TOL*maximum(sum(values(sampler.E))) : TOL
		@debug println("TOL is $(TOL), realTOL to solve for is is $(realTOL)")

		# estimate initial bias and splitting
		B = maximum(compute_bias(sampler,boundary))
    !sampler.showInfo || print_with_color(:green, sampler.ioStream, @sprintf("*** initial bias = %0.6e \n", B ) )
		splitting = ( B < realTOL/2 ) ? 1 - B/realTOL : 0.5 # update splitting parameter \in [0.5,1)
    !sampler.showInfo || print_with_color(:green, sampler.ioStream, @sprintf("*** splitting = %0.6e \n", splitting ) )

    # calculate optimal number of samples
    mySum = zeros(T,sampler.Z)
    for index::Index{d,Vector{N}} in newindexset
			mySum += ( ( sampler.Vf[index].*(sampler.W[index]*ones(sampler.Z)).^(2*λ) ).^(1/(2*λ+1) ) )
    end
    for index::Index{d,Vector{N}} in newindexset
			Nopt = ( ( sqrt(2)*erfcinv(failProb)/(splitting*realTOL) )^2 * 1/q * (sampler.Vf[index]./sampler.W[index]).^(2*λ/(2*λ+1)) .* mySum ).^(1/(2*λ))
      S[index] = ceil(N,max(3.,maximum(Nopt)))
    end
		@debug println("constant C is $(sqrt(2)*erfcinv(failProb))")
		@debug begin
  		str  = "-------------------------------------------------------------------------------- \n"
  		itype = ndims(sampler.indexSet) == 1 ? "level" : "index"
  		str *= "  "*itype*"       optimal number of samples \n"
  		str *= "-------------------------------------------------------------------------------- \n"
  		for index in sort(Set(collect(keys(S))))
    		str *= ("  $(index.indices)            "[1:13])
    		str *= @sprintf("    %d               \n",S[index])
  		end
			print(str)
		end

    # take additional samples at each level
		@debug println("begin to take samples as shown above") 
		for index::Index{d,Vector{N}} in newindexset
			samplesToTake = max( 0, S[index]-size(get(sampler.samples,index,zeros(0,0)),dir) )
      if samplesToTake > 0
        sample(sampler, samplesToTake, index )
      end
    end
		@debug inspect(sampler)
    
		# safety
		B = maximum(compute_bias(sampler,boundary))
		splitting = ( B < realTOL/2 ) ? 1 - B/realTOL : 0.5 # update splitting parameter \in [0.5,1)
		A = maximum(compute_stochastic_error(sampler, failProb, newindexset))
    @debug print_with_color(:green, sampler.ioStream, @sprintf("*** bias = %0.6e \n", B ) )
    @debug print_with_color(:green, sampler.ioStream, @sprintf("*** splitting = %0.6e \n", splitting ) )
    @debug print_with_color(:green, sampler.ioStream, @sprintf("*** stochastic error is = %0.6e \n", A ) )	
		@debug print_with_color(:red, "safety is on... ")
		if sampler.safety
			@debug print_with_color(:red, "yes\n")
      while A > splitting*realTOL
				@debug println(@sprintf("stochastic error (%0.6e) > splitting*realTOL (%0.6e), entering while loop...", A, splitting*realTOL))
        # double number of samples on index with max ratio Vest/W
        maxratio = zero(T)
        maxindex = Index(zeros(N,d))::Index{d,Vector{N}}
        for index::Index{d,Vector{N}} in newindexset
          ratio = maximum(sampler.Vest[index])/sampler.W[index]
          if ratio > maxratio
            maxratio = ratio
            maxindex = index::Index{d,Vector{N}}
          end
        end
				@debug println("the index with maximum ratio Vest/W is $(maxindex.indices)")
        n = size(sampler.samples[maxindex],dir)
				@debug println("currently $(n) samples have been taken at index $(maxindex.indices), taking an additional $(nextpow2(n+1)-n)")
        sample(sampler, nextpow2(n+1)-n, maxindex ) # round to nearest power of two
				@debug inspect(sampler)
    
				# reevaluate optimal number of samples
				# TODO maybe refactor into a function, since it is called twice
    		@debug println("reevaluating optimal number of samples...")
				for index::Index{d,Vector{N}} in newindexset
					mySum += ( ( sampler.Vf[index].*(sampler.W[index]*ones(sampler.Z)).^(2*λ) ).^(1/(2*λ+1) ) )
    		end
    		for index::Index{d,Vector{N}} in newindexset
					Nopt = ( ( sqrt(2)*erfcinv(failProb)/(splitting*realTOL) )^2 * 1/q * (sampler.Vf[index]./sampler.W[index]).^(2*λ/(2*λ+1)) .* mySum ).^(1/(2*λ))
      		S[index] = ceil(N,max(3.,maximum(Nopt)))
    		end
				@debug println("constant C is $(sqrt(2)*erfcinv(failProb))")
				@debug begin
  				str  = "-------------------------------------------------------------------------------- \n"
  				itype = ndims(sampler.indexSet) == 1 ? "level" : "index"
  				str *= "  "*itype*"       optimal number of samples \n"
  				str *= "-------------------------------------------------------------------------------- \n"
  				for index in sort(Set(collect(keys(S))))
    				str *= ("  $(index.indices)            "[1:13])
    				str *= @sprintf("    %d              \n ",S[index])
  				end
					println(str)
				end
    		
				# take these additional samples at each level
				# TODO same as above, maybe merge into a single loop
				for index::Index{d,Vector{N}} in newindexset
					samplesToTake = max( 0, S[index]-size(get(sampler.samples,index,zeros(0,0)),dir) )
      		if samplesToTake > 0
        		sample(sampler, samplesToTake, index )
      		end
    		end
				@debug inspect(sampler)
			
				B = maximum(compute_bias(sampler,boundary))
    		splitting = ( B < realTOL/2 ) ? 1 - B/realTOL : 0.5
				A = maximum(compute_stochastic_error(sampler, failProb, newindexset))	
    		@debug print_with_color(:green, sampler.ioStream, @sprintf("*** bias = %0.6e \n", B ) )
    		@debug print_with_color(:green, sampler.ioStream, @sprintf("*** splitting = %0.6e \n", splitting ) )
    		@debug print_with_color(:green, sampler.ioStream, @sprintf("*** stochastic error is = %0.6e \n", A ) )	
				@debug "I reached the end of another iteration."		
			end
			@debug print_with_color(:red, "finished safety loop!\n")
		else
			@debug print_with_color(:red, "no\n")
		end
		
    # check for convergence
		if ( L > 1 )
			@debug print("checking for convergence... ")
			error = sqrt(A^2+B^2)::T  
			converged = ( error < realTOL )
			@debug converged ? println("yes") : println("no")
			@debug println("the requested error is $(realTOL)")
			@debug println("we have error $(error)")
      # print some info to the screen
      !sampler.showInfo || print_with_color(:magenta, sampler.ioStream, 
				@sprintf("*** estimate for the bias is %0.6e \n", B ) )
      !sampler.showInfo || print_with_color(:magenta, sampler.ioStream, 
				@sprintf("*** estimate for the stochastic error is %0.6e \n", A ) )
			!sampler.showInfo || print_with_color(:magenta, sampler.ioStream, 
				@sprintf("*** estimate for the total error is %0.6e \n", error ) )
			μ = mean(sampler, newindexset)
			σ = std(sampler, newindexset)
      !sampler.showInfo || ( !converged || print_with_color(:magenta, sampler.ioStream, 
        @sprintf("*** result is %0.6e ±(%0.6e, %0.6e, %0.6e) \n", μ, σ, 2*σ, 3*σ ) ) )
		else
			@debug println("L < 1, not checking convergence... ")
		end

		# when debugging, update Lcollection
		@debug Lcollection[L+1,:] = [B, A, sqrt(A^2+B^2)]
		@debug writedlm("data/indexset_errors.txt",Lcollection)

    # print warning if no convergence and maximum level reached
    if !converged
      if L == sampler.maxL
        warn("maximum level reached and no convergence yet, sorry! :(\n")
        converged = true
      end
    end

		# if burn-in of adaptive method has passed
		if typeof(sampler.indexSet) <: AD && L == 2
			for index in boundary
				active[index]=0.
			end
			for index in oldindexset
				push!(old,index)
			end
			L = length(newindexset)
		end
   
		####
		# @debug print_with_color(:red,"*********************\n")
		# @debug print_with_color(:red,"*      WARNING      *\n")
		# @debug print_with_color(:red,"*********************\n")
		# @debug println("I set converged to false to test regression model")
	  # converged = L == sampler.maxL ? true : false	
		###

		# update L
		@debug converged || println("updating L form $(L) to $(L+1)")
		L += 1
    oldindexset = newindexset

	end
end
