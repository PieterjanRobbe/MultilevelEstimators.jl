## analyse.jl : usefull tool in analysis of MLMC use cases

# use this file to inspect behavior of expected value, variance and cost across the levels or indices 
# analyse will store the results in the data folder
# E, V, W are mean, sample variance and run time
# dE, dV, dW are 
# every column contains a header and the result for a specific direction for all levels
# directions are ordered in binary logic
function analyse(sampler;max_level=4,min_samples=20,max_samples=10000,folder=".")

	sampler.storeSamples0 = true

	d = ndims(sampler)

	# check if results are already present
	print("checking if result directories are already present... ")
	if isfile(folder*"/W.txt")
		println("yes")
		# load data
		print("reading data from files... ")
		E = readdlm(folder*"/E.txt")
		V = readdlm(folder*"/V.txt")
		W = readdlm(folder*"/W.txt")
		dE = readdlm(folder*"/dE.txt")
		dV = readdlm(folder*"/dV.txt")
		println("done")
	else
		println("no")
		# preallocate result arrays
		print("preallocating arrays... ")
		E = zeros(max_level,2^d)
		V = zeros(max_level,2^d)
		W = zeros(max_level,2^d)
		dE = zeros(max_level,2^d)
		dV = zeros(max_level,2^d)
		println("done")
	end 

	# fill index 0
	index = Index(zeros(Int64,d))
	print("checking if index $(index.indices) is empty... ")
	if W[1,1] == 0
		println("yes")
		# take samples
		nb_of_samples = max_samples
		W[1,1] = @elapsed MultilevelEstimators.sample(sampler, nb_of_samples, index)
		W[1,1] /= nb_of_samples # normalize for number of samples taken 
		# update results
		E[1,1] = squeeze(mean(sampler.samples0[index],(1,2)),(1,2))[1]
		V[1,1] = squeeze(mean(var(sampler.samples0[index],2),1),(1,2))[1]
		dE[1,1] = squeeze(mean(sampler.samples[index],(1,2)),(1,2))[1]
		dV[1,1] = squeeze(mean(var(sampler.samples[index],2),1),(1,2))[1]
		# store results
		print("saving into $(folder)/... ")
		writedlm(folder*"/E.txt",E)
		writedlm(folder*"/V.txt",V)
		writedlm(folder*"/W.txt",W)
		writedlm(folder*"/dE.txt",dE)
		writedlm(folder*"/dV.txt",dV)
		println("done")
	else
		println("no")
		println("      --> skipping index $(index.indices) (dirty)")
	end
	# loop over all other indices
	for i = 1:max_level
		for j = 2:2^d               
			index = Index(i*digits(j-1,2,d))
			print("checking if index $(index.indices) is empty... ")
			if W[i,j] == 0
				println("yes")
				# take samples
				nb_of_samples = max(min_samples, convert(Int64,round(max_samples/2^sum(index)))) 
				W[i,j] = @elapsed MultilevelEstimators.sample(sampler, nb_of_samples, Index(index))
				W[i,j] /= nb_of_samples # normalize for number of samples taken 
				# update results
				E[i,j] = squeeze(mean(sampler.samples0[index],(1,2)),(1,2))[1]
				V[i,j] = squeeze(mean(var(sampler.samples0[index],2),1),(1,2))[1]
				dE[i,j] = squeeze(mean(sampler.samples[index],(1,2)),(1,2))[1]
				dV[i,j] = squeeze(mean(var(sampler.samples[index],2),1),(1,2))[1]
				# store results
				print("saving into data/... ")
				writedlm(folder*"/E.txt",E)
				writedlm(folder*"/V.txt",V)
				writedlm(folder*"/W.txt",W)
				writedlm(folder*"/dE.txt",dE)
				writedlm(folder*"/dV.txt",dV)
				println("done")
			else
				println("no")
				println("      --> skipping index $(index.indices) (dirty)")
			end
		end
	end
end
