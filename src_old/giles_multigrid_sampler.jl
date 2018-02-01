## giles_multigrid_sampler.jl : Giles multigrid multilevel add-ons

# main sample function for Giles multigrid multilevel monte carlo
function sample_gmg(sampler, nb_of_samples, index)

	# parallel execution
	desc="Taking $(nb_of_samples) samples at index $(index.indices) "
	progress = Progress(nb_of_samples, dt=1, desc=desc, color=:black, barlen=25) # fancy progress bar
	t = @elapsed r = pmap((i)->single_sample_gmg(index,i,sampler), progress, 1:nb_of_samples)
	#t = @elapsed r = map((i)->single_sample_gmg(index,i,sampler), 1:nb_of_samples)
	the_samples_diff = reduce((x,y)->cat(2,x,y),r)

	# samples at the finest level ell
	samples_ell_diff = reshape(the_samples_diff,(size(the_samples_diff)...,1))

	# check if entry at level ell already exists
	if haskey(sampler.samples, index) # append samples
		sampler.samples[index] = hcat(sampler.samples[index], samples_ell_diff)
	else # make new entry
		sampler.samples[index] = samples_ell_diff
	end

	# update simulation time
	sampler.T[index] = t + get(sampler.T,index,0.0) # cumulative time 

	# update all other Dicts
	update_dicts(sampler, index) # on level ell

	# store samples0 if needed
	#=if sampler.storeSamples0
		the_samples = reduce((x,y)->cat(2,x[1],y[1]),r)
		samples_ell = reshape(the_samples,(size(the_samples)...,1))
		if haskey(sampler.samples0, index) # append samples
			sampler.samples0[index] = vcat(sampler.samples0[index], samples_ell)
		else # make new entry
			sampler.samples0[index] = samples_ell
		end
	end=#
end



# compute a single multigrid multilevel sample
# function to be run in parallel
# idea is that we run the sample once on every physical grid
# and return the solution after every V-cycle
# index has dimension one as physical levels, and dimension 2 as number of V-cycles
# currently only supports multi-index MC sampling
function single_sample_gmg(index, sample_no, sampler)
	ell = index[1]
	tau = index[2]

	# get the sample_no'th point from the point set generator
	# generates s x q random numbers
	xi = getPoint(sampler.numberGenerator[Index(0,0)], sample_no) # TODO fix number generator index

	# compute sample at finest grid using multigrid
	sample_ff, sample_fc = sampler.sampleFunction(xi, index, sampler)
	sample_cf, sample_cc = ell > 0 ? sampler.sampleFunction(xi, Index(ell-1,tau), sampler) : (0,0) 

	#return sample_ff, sample_ff - sample_fc - sample_cf + sample_cc            
	return sample_ff - sample_fc - sample_cf + sample_cc            
end
