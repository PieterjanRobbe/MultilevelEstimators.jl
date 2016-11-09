# add procs
addprocs(3)

# make sure the test module can be found
push!(LOAD_PATH,".")

# include the test module
using TestModule

using Base.Test

#########################################################################
# 2D ELLIPTIC SPDE WITH NEUMANN BOUNDARY CONDITIONS, MLMC, MULTIPLE QOI
#########################################################################
function test1(TOL::AbstractFloat)
	pd = 2 # physical dimension of the problem
	d = 1 # dimension of the index set (1 is multilevel)

	λ = 1.
	σ = 1.
	ν = 0.5
	p = 1.

	s = 100
	
	myIndexSet = ML() # set-up for multilevel
	myNumberGenerator = GaussianMCgenerator(s) # Monte Carlo sampler
	myMaternKernel = MaternKernel(λ,σ,ν,p)
	myGaussianFieldSampler = KLExpansion(myMaternKernel,pd,s,m0=4,maxL=6)

	myDict = Dict(
    	"indexSet" => myIndexSet,
    	"numberGenerator" => myNumberGenerator,
    	"sampleFunction" => parametrizedPDEPointEvaluation,
    	"gaussianFieldSampler" => myGaussianFieldSampler,
    	"Z" => 9
	)

	mySampler = setup(myDict)

	t = simulate(mySampler,TOL)

	E = 0.
	for index in keys(mySampler.samples)
      E += squeeze(mean(mySampler.samples[index],(1,2)),(1,2))[5]
  	end

	@test_approx_eq_eps E 0.5 TOL
end

#########################################################################
# 2D ELLIPTIC SPDE WITH DIRICHLET BOUNDARY CONDITIONS, MIQMC (TD)
#########################################################################
function test2(TOL::AbstractFloat)
	pd = 2 # physical dimension of the problem
	d = 2 # dimension of the index set

	λ = 1.
	σ = 1.
	ν = 0.5
	p = 1.

	s = 100
	q = 16
	
	myIndexSet = TD(d) # set-up for multilevel
	myNumberGenerator = GaussianQMCgenerator(s,q) # Quasi-Monte Carlo sampler
	myMaternKernel = MaternKernel(λ,σ,ν,p)
	myGaussianFieldSampler = KLExpansion(myMaternKernel,pd,s,m0=4,maxL=6)

	myDict = Dict(
    	"indexSet" => myIndexSet,
    	"numberGenerator" => myNumberGenerator,
    	"sampleFunction" => parametrizedPDEEffectiveConductivity,
    	"gaussianFieldSampler" => myGaussianFieldSampler
	)

	mySampler = setup(myDict)

	t = simulate(mySampler,TOL)
end

#########################################################################
# 2D ELLIPTIC SPDE WITH NEUMANN BOUNDARY CONDITIONS, ACMIMC
#########################################################################
function test3(TOL::AbstractFloat)
	pd = 2 # physical dimension of the problem
	d = 2 # dimension of the index set

	λ = 1.
	σ = 1.
	ν = 0.5
	p = 1.

	s = 100
	
	myIndexSet = AD(d) # set-up for multilevel
	myNumberGenerator = GaussianMCgenerator(s) # Monte Carlo sampler
	myMaternKernel = MaternKernel(λ,σ,ν,p)
	myGaussianFieldSampler = KLExpansion(myMaternKernel,pd,s,m0=4,maxL=6)

	myDict = Dict(
    	"indexSet" => myIndexSet,
    	"numberGenerator" => myNumberGenerator,
    	"sampleFunction" => parametrizedPDEEffectiveConductivity,
    	"gaussianFieldSampler" => myGaussianFieldSampler,
    	"continuate" => true,
    	"maxL" => 10
	)

	mySampler = setup(myDict)

	t = simulate(mySampler,TOL)
end

# reset random number generator for reproducibility
@everywhere srand(2016)

# run first test
test1(1e-3)

# run second test
test2(1e-3)

# run third test
test3(1e-3)