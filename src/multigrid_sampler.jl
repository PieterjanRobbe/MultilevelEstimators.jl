## multigrid_sampler.jl : multigrid multielvel add-ons

# main sample function for multigrid multilevel monte carlo
function sample_mg(sampler, nb_of_samples, index)
	is_qmc = isa(sampler.numberGenerator,QMCgenerator) 
	@assert is_qmc # MUST be true for our multigrid sampler

	# ask state of generator (for actual qmc sampling)
	state = sampler.generatorStates
	if !haskey(state,index)
		state[index] = 0
	end

	# parallel execution
	desc="Taking $(nshifts(sampler.numberGenerator)) x $(nb_of_samples) samples at level $(index[1]) "
	progress = Progress(nb_of_samples, dt=1, desc=desc, color=:black, barlen=25) # fancy progress bar
	t = @elapsed r = pmap((i)->single_sample_mg(index,i,sampler), progress, (state[index]+1):(state[index]+nb_of_samples))
	the_samples = reduce((x,y)->cat(3,x,y),r)

	# samples at the finest level ell
	samples_ell = reshape(the_samples[end,:,:]',(size(the_samples,3,2)...,1))
	# check if entry at level ell already exists
	if haskey(sampler.samples, index) # append samples
		sampler.samples[index] = vcat(sampler.samples[index], samples_ell)
	else # make new entry
		sampler.samples[index] = samples_ell
	end
	
	# when reusing samples, append them to the previous samples on all levels
	sampler.reuse && for ell = 0:(index[1]-1)
		samples_ell = reshape(the_samples[ell+1,:,:]',(size(the_samples,3,2)...,1))
		sampler.samples[Index(ell)] = vcat(sampler.samples[Index(ell)], samples_ell)
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
end

# compute a single multigrid multilevel sample
# function to be run in parallel
function single_sample_mg(index, sample_no, sampler)
	ell = index[1]
	nb_of_shifts = nshifts(sampler.numberGenerator)

	# get the sample_no'th point from the point set generator
	# generates s x q random numbers
	xi_mat = getPoint(sampler.numberGenerator, sample_no)

	# preallocate room for the sample
	shifted_sample = zeros(ell+1,nb_of_shifts)

	# loop over all shifts
	for q in 1:nb_of_shifts
		xi = xi_mat[:,q]
		# compute sample at finest grid using multigrid
		samples = sampler.sampleFunction(xi, index, sampler)
		# compute differences and store into shifted sample
		shifted_sample[1,q] = samples[1]
		shifted_sample[2:end,q] = diff(samples)
	end

	shifted_sample
end
