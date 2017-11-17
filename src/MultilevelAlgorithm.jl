#
# MULTILEVELALGORITHM.JL
#

#	Collection of functions that implement the different algorithms, MLMC, MLQMC, MIMC, MIQMC...
# Main function is "simulate", see below.

# simulate(sampler, absTol)
simulate{T<:AbstractFloat}(sampler::Sampler, absTol::T) = simulate(sampler, [absTol])

# simulate(sampler, absTol, folder)
simulate{T<:AbstractFloat}(sampler::Sampler, absTol::T, folder::AbstractString) = simulate(sampler, [absTol], folder=folder)

# simulate(sampler, absTol, failProb)
simulate{T<:AbstractFloat}(sampler::Sampler, absTol::T, failProb::T) = simulate(sampler, [absTol], failProb=failProb)

# simulate(sampler, absTol, failProb, folder)
simulate{T<:AbstractFloat}(sampler::Sampler, absTol::T, failProb::T, folder::AbstractString) = simulate(sampler, [absTol], failProb=failProb, folder=folder)

# simulate(sampler, tol, failProb=..., is_relative=..., folder=...)
# simulate takes a vector tol as input that contains the tolerances that must be solved for
function simulate{T<:AbstractFloat,N}(sampler::Sampler, tol::Vector{T}; failProb::T=0.1, is_relative::Bool=false, folder::AbstractString=".", max_time::N=typemax(Int64))

	# checks on inputs
	if any(tol .< 0) || any(isinf.(tol))
		error("supplied tolerances must be positive!")
	elseif failProb < 0 || failProb > 1
		error("failure probability failProb must be between 0 and 1!")
		#elseif !isdir(folder)
		#error("unknown folder $(folder)!")
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
	means = zeros(tol)
	errorz = zeros(length(tol),4)


	#####
	#####
	#####
	# warm-up
	tic()
	if sampler.ml_sample_fun == sample_mg
		if sampler.reuse
			sampler.ml_sample_fun(sampler,1,Index(3))
		else
			sampler.ml_sample_fun(sampler,1,Index(0))
			sampler.ml_sample_fun(sampler,1,Index(1))
			sampler.ml_sample_fun(sampler,1,Index(2))
			sampler.ml_sample_fun(sampler,1,Index(3))
		end
	end
	inspect(sampler)
	wctime[1] = toq()
	#####
	#####
	#####
	
	# run mimc 
	for i in 1:length(tol)
		# make dir
		printdir = @sprintf("%0.4e",tol[i])
		#@debug begin  
		if !isdir(folder*"/data")
			mkdir(folder*"/data")
		end
		if !isdir(folder*"/data/"*printdir)
			mkdir(folder*"/data/"*printdir)
		end
		if !isdir(folder*"/data/"*printdir*"/indexsets")
			mkdir(folder*"/data/"*printdir*"/indexsets")
		end
		#end
		# simulate
		delta_t = @elapsed (A,B) = mimc(sampler,tol[i],is_relative,failProb,folder) 
		# compute wall clock time and standard cost
		wctime[i] +=  delta_t
		for index in keys(sampler.samples)
			stcost[i] += (length(sampler.samples[index])-8)*sampler.W[index]
		end # for
		print_with_color(:cyan,"ELAPSED IS $(delta_t)\n")
		@debug begin
			print("writing wall clock times into $(folder)/data/"*printdir*"/wctime.txt... ")
			writedlm(folder*"/data/"*printdir*"/wctime.txt",cumsum(wctime))
			println("done")
			print("writing standard cost into $(folder)/data/"*printdir*"/stcost.txt... ")
			writedlm(folder*"/data/"*printdir*"/stcost.txt",stcost)
			println("done")
			print("writing tolerances into $(folder)/data/"*printdir*"/tolerances.txt... ")
			writedlm(folder*"/data/"*printdir*"/tolerances.txt",tol)
			println("done")
			# save all samples
			#      		dir_name = folder*"/data/"*printdir*"/samples"
			#      		isdir(dir_name) ? nothing : mkdir(dir_name)
			#      		for idx in keys(sampler.samples)
			#        		dir_name2 = @sprintf("%s",idx.indices)
			#        		isdir(dir_name*"/"*dir_name2) ? nothing : mkdir(dir_name*"/"*dir_name2)
			#        		for z in 1:sampler.Z
			#        			writedlm(dir_name*"/"*dir_name2*"/Z$(z).txt",sampler.samples[idx][:,:,z])
			#        		end
			#      		end
			# save number of orignal samples
						nb_of_orig_samples = [sampler.nb_of_orig_samples[ell] for ell in sort(Set(keys(sampler.nb_of_orig_samples)))]
						print("writing nb_of_orig_samples into $(folder)/data/"*printdir*"/nb_of_orig_samples.txt... ")
						writedlm(folder*"/data/"*printdir*"/nb_of_orig_samples.txt",nb_of_orig_samples)
			println("done")
			# save tol x bias x std x total error
			errorz[i,1] = tol[i]
			errorz[i,2] = A 
			errorz[i,3] = B 
			errorz[i,4] = errorz[i,2] + errorz[i,3]
			print("writing errorz into $(folder)/data/"*printdir*"/errorz.txt... ")
			writedlm(folder*"/data/"*printdir*"/errorz.txt",errorz)
			println("done")
		end # begin
		# save means
		#		for index in keys(sampler.samples)
		#			means[i] += maximum(sampler.E[index])
		#		end
		#		print("writing means into $(folder)/data/"*printdir*"/means.txt... ")
		#		writedlm(folder*"/data/"*printdir*"/means.txt",means)
		#		println("done")
		# end save means
				if cumsum(wctime)[end] > max_time
					print_with_color(:red,"wall clock time exceeded maximum time of $(max_time) secondes, aborting...\n")
					break
				end # if
	end # for

	return cumsum(wctime),cumsum(stcost),means
end # function

# actual mimc simulation
function mimc{d,T<:AbstractFloat}(sampler::Sampler{d}, TOL::T, is_relative::Bool, failProb::T, folder::AbstractString)

	# print some info to screen
	if ( sampler.showInfo )
		@printf(sampler.ioStream,"%s","--------------------------------------------------------------------------------\n")
		@printf(sampler.ioStream,"%s","*** MultilevelEstimators.jl @$(now())\n")
		@printf(sampler.ioStream,"%s","*** Simulating $(sampler.sampleFunction)\n")
		isCont = sampler.continuate ? "" : "not "
		idxSet = ( isa(sampler.indexSet,ML) || isa(sampler.indexSet,SL) || isa(sampler.indexSet,AD) ) ? "$(sampler.indexSet)" : "$(sampler.indexSet)"[1:34]
		@printf(sampler.ioStream,"%s","*** Using a $(idxSet), $(isCont)continuating \n")
		printTOL = is_relative ? "relTOL" : "absTOL"
		@printf(sampler.ioStream,"%s","*** $printTOL = "*@sprintf("%0.3e (failure probability of %0.2f)\n",TOL,failProb))
		@printf(sampler.ioStream,"%s","--------------------------------------------------------------------------------\n")
	end # if

	# loop variables
	converged = false
	L  = 0
	printdir = @sprintf("%0.4e",TOL)

	# when debugging, keep track of parameter L, bias, stochastic error, total error
	@debug Lcollection = zeros(T,0,3)

	A=0
	B=0


	# aliases
	# TODO might be incorporated into the sampler object
	N = Int64
	λ = 1/2 #TODO QMC doesn't use this anymore, so only MC, hence 1/2; sampler.numberGenerator.λ
	q = nshifts(sampler.numberGenerator[Index(zeros(Int64,ndims(sampler)))])#nshifts(sampler.numberGenerator)
	dir = isa(sampler.numberGenerator[Index(zeros(Int64,ndims(sampler)))],MCgenerator) ? 2 : 1
	is_adaptive = (typeof(sampler.indexSet) <: AD)
	@show  is_qmc = isa(sampler.numberGenerator[Index(zeros(Int64,ndims(sampler)))],QMCgenerator)

	# variable definition
	oldindexset = Set{Index{d,Vector{N}}}()

	# variables used when adaptive
	old = Set{Index{d,Vector{N}}}() # old index set
	active = Set{Index{d,Vector{N}}}() # active set
	deleted = Set{Index{d,Vector{N}}}() # deleted set

	while !converged
		# print some info to the screen
		!sampler.showInfo || print_with_color(:blue, sampler.ioStream, "*** currently running at L = $(L)...\n")

		# empty dicts
		S = Dict{Index{d,Vector{N}},N}() # dict for optimal number of samples

		# get index set in current iteration
		if !is_adaptive
			newindexset = getIndexSet(sampler.indexSet,L)::Set{Index{d,Vector{N}}}	
		elseif is_adaptive && L <= 2
			imaxgain = Index(zeros(N,d))::Index{d,Vector{N}}
			newindexset = getIndexSet(TD(d),L)::Set{Index{d,Vector{N}}}
		else
			#if L <= 2 # assume TD index set for inial estimates
			#	newindexset = getIndexSet(TD(d),L)
			#	newindexsetboundary = getBoundary(newindexset)
			#	for index in newindexsetboundary
			#		push!(active,index)
			#	end
			#	for index in setdiff(newindexset, newindexsetboundary)
			#		push!(old,index)
			#	end
			#else
			# find index with largest gain
			maxgain = zero(T)
			i = Index(zeros(N,d))::Index{d,Vector{N}}
			for index::Index{d,Vector{N}} in active
				@debug println("index $(index.indices) has gain $(sampler.P[index])")
				gain = sampler.P[index]
				if gain > maxgain
					maxgain = gain
					i = index::Index{d,Vector{N}}
				end
			end
			imaxgain = i # store variable imaxgain for further reference
			@debug println("index with maximum profit "*@sprintf("%0.6e",sampler.P[i])*" is $(i.indices)")
			delete!(active,i) # delete from active set
			push!(old,i) # add to old set and check admissables
			for p in 1:d 
				j = copy(i)::Index{d,Vector{N}}
				j[p] += 1
				if isAdmissable(old,j) # add when admissible in old set
					push!(active,j)	
				end
			end
			@debug print("the old set is an ")
			@debug prettyprint(old)
			@debug print("the active set is an ")
			@debug prettyprint(active)
			@debug print("the deleted set is an ")
			@debug prettyprint(deleted)
			newindexset = union(deleted,old, active)::Set{Index{d,Vector{N}}}
			#end
		end

		##################
		##################
		##################
#=
		myset = Set{Index{2,Vector{Int64}}}()
		for i = 0:5
			push!(myset,Index(i,1))
			push!(myset,Index(i,0))
		end

		newindexset = typeof(sampler.indexSet) <: ML ? getIndexSet(ML(),5) : myset
=#
		##################
		##################
		##################




		# print index set
		@debug print("the new index set is an ")
		@debug prettyprint(newindexset)
		# store the index set in a file
		@debug begin
			print("writing the $(length(newindexset)) indices of L = $(L) into a file...")
			# store all indices in a matrix M
			M = zeros(T, length(newindexset), d)
			cntr = 1
			for idx in sort(newindexset)
				M[cntr,1:d] = idx.indices
				cntr += 1
			end
			writedlm(folder*"/data/"*printdir*"/indexsets/L$(L).txt",M)
			println("done")
		end

		@debug begin
			print("writing the $(length(active)) active indices of L = $(L) into a file...")
			# store all indices in a matrix M
			M = zeros(T, length(active), d)
			cntr = 1
			for idx in sort(active)
				M[cntr,1:d] = idx.indices
				cntr += 1
			end
			writedlm(folder*"/data/"*printdir*"/indexsets/La$(L).txt",M)
			println("done")
		end

		@debug begin
			print("writing the $(length(deleted)) deleted indices of L = $(L) into a file...")
			# store all indices in a matrix M
			M = zeros(T, length(deleted), d)
			cntr = 1
			for idx in sort(deleted)
				M[cntr,1:d] = idx.indices
				cntr += 1
			end
			writedlm(folder*"/data/"*printdir*"/indexsets/Ld$(L).txt",M)
			println("done")
		end

		@debug if is_adaptive && L > 2
			print("writing the max gain index of L = $(L) into a file...")
			writedlm(folder*"/data/"*printdir*"/indexsets/Lg$(L).txt",imaxgain.indices)
			println("done")
		end

		# get boundary in current iteration
		boundary = Set{Index{d,Vector{N}}}()::Set{Index{d,Vector{N}}} # empty set
		if typeof(sampler.indexSet) <: ML
			push!(boundary,Index(L*ones(N,d)))
		elseif !is_adaptive
			boundary = union(boundary,getBoundary(getIndexSet(sampler.indexSet, L))::Set{Index{d,Vector{N}}})
		elseif is_adaptive && L < 3
			boundary = union(boundary,getBoundary(getIndexSet(TD(d), L))::Set{Index{d,Vector{N}}})
		else
			#boundary = union(active,deleted)
			boundary = union(active)
		end
		@debug print("the boundary is an ")
		@debug prettyprint(boundary)

		# new indices to add
		indicesToAdd = setdiff( newindexset, oldindexset )
		@debug print("the new indices are an ")
		@debug prettyprint(indicesToAdd)

		# take initial number of samples at each new index when needed...
		if L < 3000 ############# ALWAYS
			@debug println("now taking warm-up samples... ")
			for index::Index{d,Vector{N}} in indicesToAdd
				if !haskey(sampler.samples,index) # when running repeatedly, might already have samples available
					@show sampler.ml_sample_fun
					sampler.ml_sample_fun(sampler, sampler.Nstar::N, index::Index{1,Vector{N}})
				end
			end
			@debug println("done taking warm-up samples!")
			# ... or else do regression over the new indices to get estimates for variances
		else
			@debug println("now doing regression... ")
			do_regression(sampler,setdiff(indicesToAdd,Set(collect(keys(sampler.samples)))))
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
		mySum = 0.
		for index::Index{d,Vector{N}} in newindexset
			#mySum += ( ( sampler.Vf[index].*(sampler.W[index]*ones(sampler.Z)).^(2*λ) ).^(1/(2*λ+1) ) )
			mySum += maximum( sqrt.( sampler.Vf[index]*sampler.W[index] ) )
		end
		for index::Index{d,Vector{N}} in newindexset
			#Nopt = ( ( sqrt(2)*erfcinv(failProb)/(splitting*realTOL) )^2 * 1/q * (sampler.Vf[index]./sampler.W[index]).^(2*λ/(2*λ+1)) .* mySum ).^(1/(2*λ))
			Nopt = 2/(realTOL)^2 * maximum( sqrt.(sampler.Vf[index]./sampler.W[index]) ) * mySum 
			S[index] = is_qmc ? sampler.Nstar : ceil(N,max(3.,maximum(Nopt)))
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
				sampler.ml_sample_fun(sampler, samplesToTake, index )
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
				#=if sampler.ml_sample_fun == sample_mg
					C = get_covariance_matrix(sampler)
					(i,j) = ind2sub(size(C),indmax(C))
					maxindex = Index(max(i,j)-1)
					maxratio = C[i,j]
				else
					=#
					# double number of samples on index with max ratio Vest/W
					maxratio = zero(T)
					maxindex = Index(zeros(N,d))::Index{d,Vector{N}}
					for index::Index{d,Vector{N}} in newindexset
						ratio = maximum(sampler.Vest[index])/sampler.W[index]
						@debug println("index $(index.indices) has ratio $(ratio)")
						if ratio > maxratio
							maxratio = ratio
							maxindex = index::Index{d,Vector{N}}
						end
					end
					#end
					@debug println("the index with maximum ratio Vest/W is $(maxindex.indices)")
					n = size(sampler.samples[maxindex],dir)
					@debug println("currently $(n) samples have been taken at index $(maxindex.indices), taking an additional $(nextpow2(n+1)-n)")
					sampler.ml_sample_fun(sampler, nextpow2(n+1)-n, maxindex ) # round to nearest power of two
					@debug inspect(sampler)

					# reevaluate optimal number of samples
					# TODO maybe refactor into a function, since it is called twice
					@debug println("reevaluating optimal number of samples...")
					for index::Index{d,Vector{N}} in newindexset
						mySum += ( ( sampler.Vf[index].*(sampler.W[index]*ones(sampler.Z)).^(2*λ) ).^(1/(2*λ+1) ) )
					end
					for index::Index{d,Vector{N}} in newindexset
						Nopt = ( ( sqrt(2)*erfcinv(failProb)/(splitting*realTOL) )^2 * 1/q * (sampler.Vf[index]./sampler.W[index]).^(2*λ/(2*λ+1)) .* mySum ).^(1/(2*λ))
						S[index] = is_qmc ? sampler.Nstar : ceil(N,max(3.,maximum(Nopt)))
					end
					@debug println("constant C is $(sqrt(2)*erfcinv(failProb))")
					@debug begin
						str  = "-------------------------------------------------------------------------------- \n"
						itype = ndims(sampler.indexSet) == 1 ? "level" : "index"
						str *= "  "*itype*"       optimal number of samples \n"
						str *= "-------------------------------------------------------------------------------- \n"
						for index in sort(Set(collect(keys(S))))
							str *= ("  $(index.indices)            "[1:13])
							str *= @sprintf("    %d              \n",S[index])
						end
						print(str)
					end

					# take these additional samples at each level
					# TODO same as above, maybe merge into a single loop
					for index::Index{d,Vector{N}} in newindexset
						samplesToTake = max( 0, S[index]-size(get(sampler.samples,index,zeros(0,0)),dir) )
						if samplesToTake > 0
							sampler.ml_sample_fun(sampler, samplesToTake, index )
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
			if ( !is_adaptive && L > 1 ) || ( is_adaptive && L > 2 )
				@debug print("checking for convergence... ")
				error = sqrt(A^2+B^2)::T  
				#converged = ( error < realTOL )
				print(">> WARNING << using only bias to compute convergence")
				converged = B^2 < realTOL^2/2
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
				@debug println("not checking convergence... ")
			end

			###########################3
			###########################3
			###########################3

		#converged = true	
		@show converged = L == sampler.maxL
			
			###########################3
			###########################3
			###########################3
			
			# when debugging, update Lcollection
			@debug Lcollection = vcat(Lcollection,[B A sqrt(A^2+B^2)])
			@debug writedlm(folder*"/data/"*printdir*"/indexset_errors.txt",Lcollection)

			# print warning if no convergence and maximum level reached
			if !is_adaptive
				@show sampler.maxL
				@show union(newindexset,sampler.max_indexset)
				@show max_reached = any(sampler.maxL .== L) #|| !isempty(union(newindexset,sampler.max_indexset))
				if !converged && max_reached
					warn("maximum level reached and no convergence yet, sorry! :(\n")
					converged = true
				end
			else
				max_reached = false
				for index in sort(active)
					@debug println("index $(index.indices) is in max_indexset? $(in(index,sampler.max_indexset))")
					@debug println("index $(index.indices) has maxL? $(any( index.indices .== sampler.maxL))")
					if in(index,sampler.max_indexset) || any( index.indices .== sampler.maxL )
						# we should remove the index from the active set
						delete!(active,index)
						push!(old,index)
						push!(deleted,index) # push into deleted, because it still counts for the "boundary"
						warn("the index $(index.indices) exceeded maximum allowed level, removing from active set...")
					end
				end
				if length(deleted) == length(sampler.max_indexset)
					# fatal...
					converged = true 
					max_reached = true
					warn("all possible indices used in the index set, still no convergence :(")
				end
			end

			# if burn-in of adaptive method has passed
			if typeof(sampler.indexSet) <: AD && L == 2
				@debug println("adaptive and L = 2, setting up active and old sets")
				for index in setdiff(getIndexSet(TD(d),3),newindexset)
					push!(active,index)
				end
				@debug print("active set is an ")
				@debug prettyprint(active)
				for index in newindexset
					push!(old,index)
				end
				@debug print("old set is an ")
				@debug prettyprint(old)
				# regress on the active indices
				do_regression(sampler,setdiff(active,Set(collect(keys(sampler.samples)))))
			end

			# update L
			@debug converged || println("updating L form $(L) to $(L+1)")
			L += 1
			if L > 100
				error("probably something went wrong, giving up now...")
			end
			oldindexset = newindexset

		end
		return (A,B)
	end
