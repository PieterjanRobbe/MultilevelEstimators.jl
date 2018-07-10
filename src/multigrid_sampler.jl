## multigrid_sampler.jl : multigrid multielvel add-ons

# main sample function for multigrid multilevel monte carlo
function sample_mg(sampler::Sampler, nb_of_samples::N, index::Index{1,Vector{N}}) where {N}
	println("go!")
	is_qmc = isa(sampler.numberGenerator[Index(zeros(Int64,ndims(sampler)))],QMCgenerator) 
	@assert is_qmc # MUST be true for our multigrid sampler
	nb_of_shifts = nshifts(sampler.numberGenerator[Index(zeros(Int64,ndims(sampler)))])

	# ask state of generator (for actual qmc sampling)
	state = sampler.generatorStates
	if !haskey(state,index)
		state[index] = 0
	end

	# parallel execution
	desc="Taking $(nb_of_shifts) x $(nb_of_samples) samples at level $(index[1]) "
	progress = Progress(nb_of_samples, dt=1, desc=desc, color=:black, barlen=25) # fancy progress bar
	start = state[index]+1::N
	stop = state[index]+nb_of_samples::N
	#
	myf(i) = single_sample_mg(index,start:stop,i,sampler)
	myf2(i) = single_sample_mg2(index,start:stop,i,sampler)
	println(desc)
	if sampler.reuse
		t = @elapsed r = pmap(myf, 1:nb_of_shifts)
	else
		t = @elapsed r = pmap(myf2, 1:nb_of_shifts)
	end
	#
	#pmap((i)->svd(randn(i)),1:10)
	#t = @elapsed r = pmap(single_sample_mg, fill(index,stop-start+1), start:stop, fill(sampler,stop-start+1) )
	the_samples = reduce((x,y)->cat(3,x,y),r)#::Array{Float64,3}
	the_samples = permutedims(the_samples,[1,3,2])

	# samples at the finest level ell
	samples_ell = reshape(the_samples[end,:,:]',(size(the_samples,3,2)...,1))
	# check if entry at level ell already exists
	if haskey(sampler.samples, index) # append samples
		#@everywhere gc()
		sampler.samples[index] = vcat(sampler.samples[index], samples_ell)
		#@everywhere gc()
	else # make new entry
		#@everywhere gc()
		sampler.samples[index] = samples_ell
		#@everywhere gc()
	end

	if haskey(sampler.nb_of_orig_samples, index) # append samples
		sampler.nb_of_orig_samples[index] += nb_of_samples
	else
		sampler.nb_of_orig_samples[index] = nb_of_samples
	end
	
	# when reusing samples, append them to the previous samples on all levels
	sampler.reuse && for ell = 0:(index[1]-1)
		samples_ell = reshape(the_samples[ell+1,:,:]',(size(the_samples,3,2)...,1))
		if haskey(sampler.samples,Index(ell))
			sampler.samples[Index(ell)] = vcat(sampler.samples[Index(ell)], samples_ell)
		else
			sampler.samples[Index(ell)] = samples_ell
		end
	end

	# update generator state
	state[index] += nb_of_samples

	# update simulation time
	sampler.T[index] = t + get(sampler.T,index,0.0) # cumulative time 

	# update all other Dicts
	update_dicts(sampler, index) # on level ell
	sampler.reuse && for ell = 0:(index[1]-1) # on all other levels
		update_dicts(sampler, Index(ell))
	end

	return nothing
end

# compute a single multigrid multilevel sample
# function to be run in parallel
function single_sample_mg(index::Index{1,Vector{N}}, sample_nos, shift_no::N, sampler::Sampler) where {N}
	
	ell = index[1]
	nb_of_shifts = nshifts(sampler.numberGenerator[Index(zeros(Int64,ndims(sampler)))])

	# get the sample_no'th point from the point set generator
	# generates s x q random numbers
	#xi_mat = getPoint(sampler.numberGenerator[index], sample_no)
	#xi_mat = randn(ndims(sampler.numberGenerator[Index(0)]),nb_of_samples)

	# preallocate room for the sample
	shifted_sample = zeros(ell+1,length(sample_nos))#nb_of_shifts)

	# loop over all shifts
	
	for q in 1:length(sample_nos)#1:nb_of_shifts
		xi = getPoint(sampler.numberGenerator[index], sample_nos[q])[:,shift_no]
		#xi = xi_mat[:,q]::Vector{Float64}
		# compute sample at finest grid using multigrid
		samples = sampler.sampleFunction(xi, index, sampler)
		# compute differences and store into shifted sample
		shifted_sample[1,q] = samples[1]
		if ell > 0
			shifted_sample[2:end,q] = diff(samples)
		end
	end
	
	shifted_sample::Matrix{Float64}
end

function single_sample_mg2(index::Index{1,Vector{N}}, sample_nos, shift_no::N, sampler::Sampler) where {N}
	
	ell = index[1]
	nb_of_shifts = nshifts(sampler.numberGenerator[Index(zeros(Int64,ndims(sampler)))])

	# preallocate room for the sample
	shifted_sample = zeros(1,length(sample_nos))#nb_of_shifts)

	# loop over all shifts
	

	for q in 1:length(sample_nos)#1:nb_of_shifts
		xi = getPoint(sampler.numberGenerator[index], sample_nos[q])[:,shift_no]
		#xi = xi_mat[:,q]::Vector{Float64}
		# compute sample at finest grid using multigrid
		samples = sampler.sampleFunction(xi, index, sampler)
		# compute differences and store into shifted sample
		shifted_sample[1,q] = samples[1]
		if ell > 0
			samples = sampler.sampleFunction(xi, Index(index[1]-1), sampler)
			shifted_sample[1,q] -= samples
		end
	end
	
	shifted_sample::Matrix{Float64}
end
