## prashant_algorithm.jl : implementation of Prshant's MGMLMC algorithm

# continuation sampling
function prashant_simulate(sampler, tol, folder)
	# tolerances to solve for
	tols = 1.5.^(sampler.nTOL-(1:sampler.nTOL))*tol

	# preallocate values that will be collected
	means = zeros(length(tols))
	runtime = zeros(length(tols))

	# run repeated MIMC simulation
	L = 0
	for t in 1:length(tols)
		## run and time
		delta_t = @elapsed prashant_mgmlmc(sampler, tols[t], L)
		runtime[t] = delta_t
		## make dir structure
		printdir = @sprintf("%0.4e",tols[t])
		if !isdir(folder*"/data")
			mkdir(folder*"/data")
		end
		if !isdir(folder*"/data/"*printdir)
			mkdir(folder*"/data/"*printdir)
		end
		## collect mean values
		for index in keys(sampler.samples)
			means[t] += maximum(sampler.E[index])
		end
		L = length(sampler.samples)-1
		## write means
		print("writing means into $(folder)/data/"*printdir*"/means.txt... ")
		writedlm(folder*"/data/"*printdir*"/means.txt",means)
		println("done")
		## write runtime
		print("writing times into $(folder)/data/"*printdir*"/wctime.txt... ")
		writedlm(folder*"/data/"*printdir*"/wctime.txt",cumsum(runtime))
		println("done")
		# save number of orignal samples
		nb_of_orig_samples = [sampler.nb_of_orig_samples[ell] for ell in sort(Set(keys(sampler.nb_of_orig_samples)))]
		print("writing nb_of_orig_samples into $(folder)/data/"*printdir*"/nb_of_orig_samples.txt... ")
		writedlm(folder*"/data/"*printdir*"/nb_of_orig_samples.txt",nb_of_orig_samples)
		println("done")
		print_with_color(:red,@sprintf("ELAPSED = %0.2f\n",cumsum(runtime)[t]))
	end
	return means
end

# Prashant's algorithm
function prashant_mgmlmc(sampler, tol, L=1)
	#L = 1
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
				prashant_sample_mg(sampler, sampler.Nstar, index)
			end
		end
		@debug inspect(sampler)

		## estimate sample variances
		# done automatically, access the variance as sampler.Vf

		## switch between prashant and cauchy schwarz for optimal number of samples
		if !sampler.use_batches
			if !sampler.is_cauchy_schwarz
				## compute optimal total number of samples
				@debug println("Computing optimal number of samples...")
				mysum = 0.
				for index in indexset
					mysum += ( maximum(sampler.Vf[index])*(sampler.W[index]) )^(1/2) 
				end
				for index in indexset
					S[index] = ceil( 2/tol^2 * ( maximum(sampler.Vf[index])/sampler.W[index] )^(1/2) * mysum )
				end
			else
				## cauchy-schwarz estimate for optimal number of samples
				@show nb_of_samples = cauchy_schwarz_opt_samples(sampler, tol)/2
				for index in indexset
					S[index] = ceil(nb_of_samples[index[1]+1])
				end
			end
		
			## compute optimal number of MG samples
			if sampler.reuse
				for index in indexset
					if !(index.indices[1] == L)
						S_bar[index] = S[index] - S[Index(index.indices[1]+1)]
					else
						S_bar[index] = S[index]
					end
				end
			else
				S_bar = S
			end
			for index in sort(indexset)
				@debug println("optimal number of samples at level $(index[1]) is    $(S[index])")
				sampler.reuse && @debug println("optimal number of MG samples at level $(index[1]) is $(S_bar[index])")
			end

			## take the additionally required samples
			@debug println("Start taking additional samples...")
			for index in indexset
				samples_to_take = max( 0, S_bar[index] - size(sampler.samples[index],2))
				samples_to_take == 0 ? nothing : prashant_sample_mg(sampler, samples_to_take, index)
			end
			@debug inspect(sampler)
		end

		## batching
		if sampler.use_batches || sampler.safety
			cov_mat = comp_cov_mat(sampler, indexset)
			var_est = sum(cov_mat)
			while var_est > tol^2/2
				cost_mat = comp_cost_mat(sampler, indexset)
				(i,j) = ind2sub((L+1,L+1),indmax(cov_mat./sqrt(cost_mat)))
				index = Index(max(i,j)-1)
				samples_already_taken = length(sampler.samples[index])
				samples_to_take = nextpow2(samples_already_taken+1)-samples_already_taken
				prashant_sample_mg(sampler, samples_to_take, index)
				cov_mat = comp_cov_mat(sampler, indexset)
				@debug show(IOContext(STDOUT, limit=true), "text/plain", cov_mat)
				@debug print("\n")
				@debug show(IOContext(STDOUT, limit=true), "text/plain", cov_mat./sqrt(cost_mat))
				@debug print("\n")
				var_est = sum(cov_mat)
				@debug inspect(sampler)	
			end
		end

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
		#=	var_est = 0.
			for index in indexset
				var_est += maximum(sampler.Vest[index])
  			end
			println("variance of estimator is $(var_est), and should be smaller then $(tol^2/2): $(var_est<tol^2/2)")
		=#end
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

## Cauchy Scwarz residual equation
function cauchy_schwarz_residu!(x, fx, V, C, ϵ)
	# extract information
	L = length(V) - 1
	N = x[1:L+1]
	λ= x[L+2]

	# compute residuals
	mysqrt = sqrt.(V./N)
	for ell = 0:L
		idx = ell+1
		fx[idx] = C[idx] + λ/N[idx]^2 * ( V[idx] - sqrt(N[idx]) * sum( mysqrt ) )
	end
	fx[L+2] = sum(mysqrt*mysqrt')-ϵ^2/2
end

## Find optimal number of samples using Cauchy-Schwarz
function cauchy_schwarz_opt_samples(sampler, epsilon)
	# extract V and C
	L = length(sampler.Vf)-1
	V = zeros(L+1)
	C = zeros(L+1)
	for ell = 0:L
		V[ell+1] = maximum(sampler.Vf[Index(ell)])
		C[ell+1] = sampler.W[Index(ell)]
	end

	# solve nonlinear system
	x0 = append!(ceil( 2/epsilon^2 * sqrt.(V./C) * sum( sqrt.( V.*C) ) ), 1e10)
	try
		@repeat 20 try
			println("check")
			res = nlsolve( (x,fx) -> cauchy_schwarz_residu!(x,fx,V,C,epsilon), x0 )
			return res.zero[1:end-1]
		catch e
			@retry if e == DomainError()
				x0[end] *= 10 # retry with bigger lambda
			end
		end
	catch e
		if !(e == DomainError() && all(x0[1:end-1].==1.0))
			rethrow(e)
		else
			return x0[1:end-1]
		end
	end
end

## Compute covariance matrix using batching algorithm
function comp_cov_mat(sampler, indexset)
	L = length(indexset) - 1
	C = zeros(L+1,L+1)
	for i in 1:L+1
		ell = Index(i-1)
		for j = 1:L+1
			tau = Index(j-1)
			# determine number of batches
			n1 = length(sampler.samples[ell])
			n2 = length(sampler.samples[tau])
			nbatches = floor(Int64,sqrt(min(n1,n2)))
			batchsize1 = floor(Int64,n1/nbatches)
			batchsize2 = floor(Int64,n2/nbatches)
			# compute means
			x = [mean(sampler.samples[ell][(b-1)*batchsize1+1:b*batchsize1]) for b = 1:nbatches-1]
			push!(x,mean(sampler.samples[ell][(nbatches-1)*batchsize1+1:end]))
			y = [mean(sampler.samples[tau][(b-1)*batchsize2+1:b*batchsize2]) for b = 1:nbatches-1]
			push!(y,mean(sampler.samples[tau][(nbatches-1)*batchsize2+1:end]))
			# compute covariance
			C[i,j] = cov(x,y)
		end
	end
	return C
end

function comp_cost_mat(sampler, indexset)
	W = []
	for ell in sort(indexset)
		push!(W,sampler.W[ell])
	end
	return W*W'
end

