## prashant_algorithm.jl : implementation of Prshant's MGMLMC algorithm

# continuation sampling
function prashant_simulate(sampler, tol, folder)
	# tolerances to solve for
	tols = 1.5.^(sampler.nTOL-(1:sampler.nTOL))*tol

	# preallocate values that will be collected
	means = zeros(length(tols))
	runtime = 0.0

	# run repeated MIMC simulation
	for t in 1:length(tols)
		# run and time
		delta_t = @elapsed prashant_mgmlmc(sampler, tols[t])
		runtime += delta_t
		# collect and save mean value
		for index in keys(sampler.samples)
			means[t] += maximum(sampler.E[index])
		end
		printdir = @sprintf("%0.4e",tols[t])
		if !isdir(folder*"/data")
			mkdir(folder*"/data")
		end
		if !isdir(folder*"/data/"*printdir)
			mkdir(folder*"/data/"*printdir)
		end
		print("writing means into $(folder)/data/"*printdir*"/means.txt... ")
		writedlm(folder*"/data/"*printdir*"/means.txt",means)
		println("done")
		print_with_color(:red,@sprintf("ELAPSED = %0.2f\n",runtime))
	end
	return means
end

# Prashant's algorithm
function prashant_mgmlmc(sampler, tol)
	L = 1
	converged = false

	while !converged
		@debug print_with_color(:cyan,"--------------------\n")
		@debug print_with_color(:cyan,@sprintf("        L = %i\n",L))
		@debug print_with_color(:cyan,"--------------------\n")
		
		indexset = getIndexSet(ML(),L)
		S = Dict{Index{1,Vector{Int64}},Int64}() # dict for optimal number of samples
		S_bar = Dict{Index{1,Vector{Int64}},Int64}() # dict for optimal number of MG samples

		## take warm-up samples on all required levels
		@debug println("Start taking warm-up samples...")
		for index in indexset
			if !haskey(sampler.samples,index)
				sampler.ml_sample_fun(sampler, sampler.Nstar, index)
			end
		end
		@debug inspect(sampler)

		## estimate sample variances
		# done automatically, access the variance as sampler.Vf

		## compute optimal total number of samples
		mysum = 0.
		for index in indexset
			mysum += ( maximum(sampler.Vf[index])*(sampler.W[index]) )^(1/2) 
		end
		for index in indexset
			S[index] = ceil( 2/tol^2 * ( maximum(sampler.Vf[index])/sampler.W[index] )^(1/2) * mysum )
		end

		## compute optimal number of MG samples
		@debug println("Computing optimal number of samples...")
		for index in indexset
			if !(index.indices[1] == L)
				S_bar[index] = S[index] - S[Index(index.indices[1]+1)]
			else
				S_bar[index] = S[index]
			end
		end
		for index in sort(indexset)
			@debug println("optimal number of samples at level $(index[1]) is    $(S[index])")
			@debug println("optimal number of MG samples at level $(index[1]) is $(S_bar[index])")
		end

		## take the additionally required samples
		@debug println("Start taking additional samples...")
		for index in indexset
			samples_to_take = max( 0, S_bar[index] - size(sampler.samples[index],2))
			samples_to_take == 0 ? nothing : sampler.ml_sample_fun(sampler, samples_to_take, index)
		end
		@debug inspect(sampler)

		## test for convergence
		if L > 1
			@debug println("Testing for convergence")
			(A, alpha) = leastSquaresFit(sampler.E, 1)
			@debug println("αis $(alpha[1])")
			converged = maximum(sampler.E[Index(L)]) < (2^alpha[1]+1)*tol/sqrt(2)
			@debug println("$(maximum(sampler.E[Index(L)])) < $((2^alpha[1]+1)*tol/sqrt(2))? $(converged ? "yes" : "no")")
		end
		if !converged
			L += 1
		end
		@debug converged && begin
			print_with_color(:blue,"I converged!\n")
			var_est = 0.
			for index in indexset
				var_est += maximum(sampler.Vest[index])
  			end
			println("variance of estimator is $(var_est), and should be smaller then $(tol^2/2): $(var_est<tol^2/2)")
		end
	end
end

# main sample function for multigrid multilevel monte carlo
function prashant_sample_mg(sampler, nb_of_samples, index)

	# block distribution
	BLOCK_SIZE = sampler.Nstar
	nblocks = convert(Int64,ceil(nb_of_samples/BLOCK_SIZE))

	# parallel execution
	desc="Taking $(nb_of_samples) samples at level $(index[1]) "
	progress = Progress(nblocks, dt=1, desc=desc, color=:black, barlen=25) # fancy progress bar
	t = @elapsed r = pmap((block)->prashant_single_sample_mg(index,sampler,BLOCK_SIZE), progress, 1:nblocks)
	the_samples = reduce(vcat,r)
	

	# samples at the finest level ell
	samples_ell = reshape(the_samples[1:nb_of_samples,end],(1,nb_of_samples,1))
	# check if entry at level ell already exists
	if haskey(sampler.samples, index) # append samples
		sampler.samples[index] = hcat(sampler.samples[index], samples_ell)
		sampler.nb_of_orig_samples[index] += nb_of_samples
	else # make new entry
		sampler.samples[index] = samples_ell
		sampler.nb_of_orig_samples[index] = nb_of_samples
	end
	
	# when reusing samples, append them to the previous samples on all levels
	sampler.reuse && for ell = 0:(index[1]-1)
		samples_ell = reshape(the_samples[1:nb_of_samples,ell+1],(1,nb_of_samples,1))
		sampler.samples[Index(ell)] = hcat(sampler.samples[Index(ell)], samples_ell)
	end

	# update simulation time
	sampler.T[index] = t + get(sampler.T,index,0.0) # cumulative time 

	# update all other Dicts
	update_dicts(sampler, index) # on level ell
	sampler.reuse && for ell = 0:(index[1]-1) # on all other levels
		update_dicts(sampler, Index(ell))
	end
end

# compute a single multigrid multilevel sample
# function to be run in parallel
function prashant_single_sample_mg(index, sampler,BLOCK_SIZE)
	ell = index[1]

	samples = zeros(BLOCK_SIZE, index[1]+1)
	for i = 1:BLOCK_SIZE

		# generates s random numbers
		xi = randn(ndims(sampler.numberGenerator[Index(0)]))

		# compute sample at finest grid using multigrid
		samples[i,:] = sampler.sampleFunction(xi, index, sampler)
	
		# compute differences and store into shifted sample
		samples[i,2:end] = diff(samples[i,:])

	end

	return samples
end
