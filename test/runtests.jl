# add procs
addprocs(3)

# set RNG for reproducability
@everywhere srand(2016)

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
	p = 1
	rho(x,y) = matern(λ,σ,ν,p,x,y)

	s = 100
	
	indexset = createIndexSet(ML,d)
	numberGenerator = GaussianMCgenerator(d,s)
	gaussianFieldSampler = createKLexpansion(pd,λ,σ,ν,s,cov=rho)
	sampleFunction = parametrizedPDEpointEvaluation
	
	mySettings = Settings(indexset, numberGenerator, gaussianFieldSampler, sampleFunction, Z=9)

	mySampler = createSampler(d,mySettings)

	(E, V, A, B, splitting, S) = simulate(mySampler,TOL)

	@test_approx_eq_eps E[5] 0.5 TOL

	@test maximum(A + B) <= TOL
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
	p = 1
	rho(x,y) = matern(λ,σ,ν,p,x,y)

	s = 100
	
	q = 16

	indexset = createIndexSet(AD,d)
	numberGenerator = GaussianQMCgenerator(d,s,q)
	gaussianFieldSampler = createKLexpansion(pd,λ,σ,ν,s,cov=rho)
	sampleFunction = parametrizedPDEEffectiveConductivity

	mySettings = Settings(indexset, numberGenerator, gaussianFieldSampler, sampleFunction)

	mySampler = createSampler(d,mySettings)

	(E, V, A, B, splitting, S) = simulate(mySampler,TOL)

	@test maximum(A + B) <= TOL
end

# run first test
test1(1e-1)
@time test1(1e-3)

# run second test
test1(1e-1)
@time test2(1e-3)