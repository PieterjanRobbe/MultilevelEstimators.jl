var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#MultilevelEstimators.jl-Documentation-1",
    "page": "Home",
    "title": "MultilevelEstimators.jl Documentation",
    "category": "section",
    "text": "For installation instructions, see Installation.For an example on how to use this package, see Example.A full description of the functionality is described in the Manual."
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "MultilevelEstimators is not added to the Julia package manager (just yet), but, the package can easily be installed by cloning the git repository. From the Julia REPL, type ] to enter Pkg mode and runpkg> add https://github.com/PieterjanRobbe/MultilevelEstimators.jl.gitThis will install the main functionality.note: Note\nThe following packages are optional.For automatic generation of reports and figures you will need the Reporter.jl  packagepkg> add https://github.com/PieterjanRobbe/Reporter.jl.gitFinally, to run the example problems, you can install the Multigrid solvers and the package to solve lognormal diffusion problems:pkg> add https://github.com/PieterjanRobbe/SimpleMultigrid.jl\n[...]\n\npkg> add https://github.com/PieterjanRobbe/NotSoSimpleMultigrid.jl\n[...]\n\npkg> add https://github.com/PieterjanRobbe/LognormalDiffusionProblems.jl\n[...]"
},

{
    "location": "#Features-1",
    "page": "Home",
    "title": "Features",
    "category": "section",
    "text": "This package featuresImplementation of Multilevel and Multi-Index (Quasi-)Monte Carlo methods.\nNative support for parallelization on multicore machines.\nFull control of advanced algorithm options such as variance regression, mean square error splitting, ...\nAutomatic generation of reports with print-quality figures using Reporter.jl.Most recent addition to the package is unbiased Multilevel and Multi-Index Monte Carlo with automatic learning of the discrete distribution of samples across all levels or mult-indices. (See index set type U.)"
},

{
    "location": "#References-1",
    "page": "Home",
    "title": "References",
    "category": "section",
    "text": "The algorithms implemented in this package are loosely based on the following papers:Basic Multilevel Monte Carlo method:Giles, M. B. Multilevel Monte Carlo Path Simulation. Operations Research 56.3 (2008): 607-617.Variance regression and continuation:Collier, N., Haji-Ali, A. L., Nobile, F., von Schwerin, E., and Tempone, R. A Continuation Multilevel Monte Carlo Algorithm. BIT Numerical Mathematics 55.2 (2015): 399-432. Quasi-Monte Carlo and Multilevel Quasi-Monte Carlo methods:Giles, M. B., and Waterhouse, B. J. Multilevel Quasi-Monte Carlo Path Simulation. Advanced Financial Modelling, Radon Series on Computational and Applied Mathematics (2009): 165-181.Multi-Index (Quasi-)Monte Carlo methods:Robbe, P., Nuyens, D., and Vandewalle, S. A Multi-Index Quasi-Monte Carlo Algorithm for Lognormal Diffusion Problems. SIAM Journal on Scientific Computing 39.5 (2017): S851-S872.Adaptive Multi-Index methods:Robbe, P., Nuyens, D., and Vandewalle, S. A Dimension-Adaptive Multi-Index Monte Carlo Method Applied to a Model of a Heat Exchanger. International Conference on Monte Carlo and Quasi-Monte Carlo Methods in Scientific Computing. Springer Proceedings in Mathematics & Statistics 241 (2018): 429-445.Unbiased estimation:Robbe, P., Nuyens, D. and Vandewalle, S. Recycling Samples in the Multigrid Multilevel (Quasi-)Monte Carlo Method. SIAM Journal on Scientific Computing, to appear (2019).\nRobbe, P., Nuyens, D. and Vandewalle, S. Enhanced Multi-Index Monte Carlo by means of Multiple Semi-coarsened Multigrid for Anisotropic Diffusion Problems. In preparation (2019)."
},

{
    "location": "example/#",
    "page": "Example",
    "title": "Example",
    "category": "page",
    "text": ""
},

{
    "location": "example/#Example-1",
    "page": "Example",
    "title": "Example",
    "category": "section",
    "text": "As a benchmark problem, we consider the elliptic PDE with random coefficients- nabla ( a(x omega) nabla u(x omega) ) = f(x)defined on the unit squre D = 0 1^2, with homogeneous boundary conditions u(x cdot) = 0 on partial D. The uncertain diffusion coefficient is given as a lognormal random field  log a(x omega) = z(x omega) where z(x omega) is a Gaussian random field (see below).note: Note\nThe following examples assume you have already installed the test dependencies as outlined in the Installation guide."
},

{
    "location": "example/#Lognormal-diffusion-problems-1",
    "page": "Example",
    "title": "Lognormal diffusion problems",
    "category": "section",
    "text": "First, we need to define a random field that will be used for the uncertain diffusion coefficient of the PDE. Fortunately, a Julia package that generates and samples from such random fields already exists.using GaussianRandomFieldsTo specify the spatial correlation in the random field, we consider the Mat\\\'ern covariance function with a given correlation length and smoothness. This function is predefined in the package. using GaussianRandomFields # hide\ncorr_length = 0.5\nsmoothness = 2.0\ncov_fun = CovarianceFunction(2, Matern(corr_length, smoothness))We generate samples of this smooth random field using a Karhunen-Lo`eve expansion with 250 terms. To this end, we discretize the domain D into a regular grid with grid size Delta x = Delta y = 1/128.n = 128\npts = range(0, stop=1, length=n)\ngrf = GaussianRandomField(cov_fun, KarhunenLoeve(250), pts, pts)A sample of this Gaussian random field is returned by either calling the sample-function directly, or by supplying the random numbers used for sampling from the random field.sample(grf)\nsample(grf, xi = randn(250))The keyword xi corresponds to the sample omega in the defintion of the PDE.We vizualize a sample of this random field.using Plots\ncontourf(sample(grf))<img src=\"../assets/grf.png\" width=\"400\">Many more options and details regarding the generation of Gaussian random fields can be found in the tutorial that accomagnies the package GaussianRandomFields.jl. Suppose we are interested in computing the expected value of a quantity of interest derived from the solution of the PDE. For example, we might be interested in the value of the solution u(x omega) at the point x=y=05."
},

{
    "location": "example/#Multilevel-Monte-Carlo-1",
    "page": "Example",
    "title": "Multilevel Monte Carlo",
    "category": "section",
    "text": "The basics of any multilevel method is a hierarchy of approximations of the model problem with an increasing accuracy, but corresponding increasing cost. For the lognormal diffusion example, such a hierarchy is provided by solving the PDE on an ever finer grid. We define a hierarchy of Gaussian random field generators usinggrfs = [GaussianRandomField(C, KarhunenLoeve(250), 1/i:1/i:1-1/i, 1/i:1/i:1-1/i) for i in 2 .^(2:8)]The corresponding system matrix that results from a finite difference discretization of the PDE is returned by the elliptic2d-function from the SimpleMultigrid package.z = sample(grf)\na = exp(z)\nA = elliptic2d(a)"
},

{
    "location": "example/#Creating-an-Estimator-1",
    "page": "Example",
    "title": "Creating an Estimator",
    "category": "section",
    "text": "The main type provided by MultilevelEstimators.jl is an Estimator. This type has two parametric subtypes, Estimator{<:AbstractIndexSet, <:AbstractSampleMethod}, where AbstractIndexSet is an abstract type for the  index set, and AbstractSampleMethod is an abstract type for the sample method. For example, a Multilevel Monte Carlo method will have type signature Estimator{ML, MC}. For a complete list of possible index set types, see IndexSet. For a complete list of possible sample method types, see SampleMethod.The main user interaction required by MultilevelEstimators.jl is the definition of a sample function. This function must return a sample of the quantity of interest, and of its difference, for the given discretization parameter (or level) and random parameters. Here is an example for the lognormal diffusion problem:function sample_lognormal(level::Level, \\omega::Vector{<:Real}, grf::GaussianRandomField)\n\n    # solve on finest grid\n    z  = sample(grf, xi=\\omega)\n    af = exp(z)\n    Af = elliptic2d(af)\n    bf = fill(eltype(Af), size(Af, 1))\n    uf = Af\\bf\n    Qf = uf[length(uf) >> 1]\n\n    if level == 0\n        Qc = zero(Qf)\n    else\n        ac = view(af, 2:2:end, 2:2:end)\n        Ac = elliptic2d(ac)\n        bc = fill(eltype(Af), size(Ac, 1))\n        uc = Ac\\bc\n        Qc = uc[length(uc) >> 1]\n    end\n    Qf - Qc\nendWhen the level parameter is zero, we return the value of the quantity of interest on the coarsest mesh. For a larger level parameter, we return both the difference of the quantity of interest with a coarser grid solution, and the value of the quantity of interest itself.The sample function can only have two input parameters: a level and a vector of random parameters. For sample functions that require data (such as the precomputed Gaussian random fields in the example above), one might need to add a convenience function:sample_lognormal(level, \\omega) = sample_lognormal(level, \\omega, grfs[level + one(level)])Finally, for the construction of an Estimator, we also need to specify the number of random parameters and their distribution. For the lognormal diffusion problem, these random parameters are normally distributed.distributions = [Normal() for i in 1:250]A full list of predefined distributions can be found in the manual for Distributions.An estimator for the lognormal diffusion problem can thus be created usingestimator = Estimator(ML(), MC(), sample_lognormal, distributions)"
},

{
    "location": "example/#Specifying-options-1",
    "page": "Example",
    "title": "Specifying options",
    "category": "section",
    "text": "Different options and settings can be passed on to the Estimator by supplying the appropriate keyword argument. For example, to set the number of warm up samples to 10, callestimator = Estimator(ML(), MC(), sample_lognormal, distributions, nb_of_warm_up_samples = 10)A complete list of optional arguments can be found in the manual for Estimator."
},

{
    "location": "example/#Running-a-simulation-1",
    "page": "Example",
    "title": "Running a simulation",
    "category": "section",
    "text": "To start a simulation that computes the expected value of the quantity of interest up to an absolute (root mean square) error of 1e-3, callh = run(estimator, 1e-3)This function will return a History object that contains usefull diagnostics about the simulation. See History for a complete list of diagnostic entries."
},

{
    "location": "example/#Vizualization-of-the-result-using-[Reporter](@ref)-1",
    "page": "Example",
    "title": "Vizualization of the result using Reporter",
    "category": "section",
    "text": "MultilevelEstimators.jl can automatically build a set of diagnostic figures based on a History file. Load the package byusing Reporterand generate a report by callingreport(h)If all goes well, you should now see a webpage with diagnostic information about the simulation, similar to this one. Print-ready .tex-figures are stored locally under figures/."
},

{
    "location": "example/#Multilevel-Quasi-Monte-Carlo-1",
    "page": "Example",
    "title": "Multilevel Quasi-Monte Carlo",
    "category": "section",
    "text": "Quasi-Monte Carlo is an alternative way to pick the samples omega. The random samples in the Monte Carlo method are replaced by deterministically well-chosen points that increase the accuracy of the estimation with respect to the number of samples. A well-known type of such point sets are `rank-1 lattice rules, implemented here as LatticeRule32. These points are used by default when calling the QMC-version of the Estimatorestimator = Estimator(ML(), QMC(), sample_lognormal, distributions)Under the hood, the default lattice rule uses a 3600-dimensional generating vector from Kuo et al.. It is easy to provide your own generating vector with the optional point_generator-argument.estimator = Estimator(ML(), QMC(), sample_lognormal, point_generator=LatticeRule32(\"my_gen_vec.txt\"))The number of shifts is controlled by the nb_of_shifts-keyword. We recommand a value between 10 (default) and 30. The number of samples is updated carefully, with a default multiplication factor of 1.2 (instead of the standard and theoretically justifiable value 2, which we found to be too aggressive). This number can be controlled with the keyword sample_mul_factor. "
},

{
    "location": "example/#Multi-Index-Monte-Carlo-1",
    "page": "Example",
    "title": "Multi-Index Monte Carlo",
    "category": "section",
    "text": ""
},

{
    "location": "example/#Multi-index-sets-1",
    "page": "Example",
    "title": "Multi-index sets",
    "category": "section",
    "text": "using MultilevelEstimators # hide\nindex_set = TD(2)\nfor L = 0:5\n    println(\"level: \", L)\n    print(index_set, L)\nendlognormal example"
},

{
    "location": "example/#Adaptive-Multi-index-Monte-Carlo-1",
    "page": "Example",
    "title": "Adaptive Multi-index Monte Carlo",
    "category": "section",
    "text": "with penalization"
},

{
    "location": "example/#Unbiased-estimation-1",
    "page": "Example",
    "title": "Unbiased estimation",
    "category": "section",
    "text": ""
},

{
    "location": "example/#Parallelization-1",
    "page": "Example",
    "title": "Parallelization",
    "category": "section",
    "text": ""
},

{
    "location": "manual/#",
    "page": "Manual",
    "title": "Manual",
    "category": "page",
    "text": ""
},

{
    "location": "manual/#Manual-1",
    "page": "Manual",
    "title": "Manual",
    "category": "section",
    "text": ""
},

{
    "location": "manual/#MultilevelEstimators.Level",
    "page": "Manual",
    "title": "MultilevelEstimators.Level",
    "category": "type",
    "text": "Level(l::Integer)\n\nReturn a level.\n\nExamples\n\njulia> level = Level(2)\n2\n\nSee also: Index\n\n\n\n\n\n"
},

{
    "location": "manual/#MultilevelEstimators.Index",
    "page": "Manual",
    "title": "MultilevelEstimators.Index",
    "category": "type",
    "text": "Index(i::Integer...)\n\nReturn a multi-index.\n\nExamples\n\njulia> index = Index(2, 1)\n(2, 1)\n\nSee also: Level\n\n\n\n\n\n"
},

{
    "location": "manual/#Index-1",
    "page": "Manual",
    "title": "Index",
    "category": "section",
    "text": "LevelIndex"
},

{
    "location": "manual/#MultilevelEstimators.SL",
    "page": "Manual",
    "title": "MultilevelEstimators.SL",
    "category": "type",
    "text": "SL()\n\nReturn a single-level index set.\n\nExamples\n\njulia> SL()\nSL\n\nSee also: ML, FT, TD, HC, ZC, AD, U\n\n\n\n\n\n"
},

{
    "location": "manual/#MultilevelEstimators.ML",
    "page": "Manual",
    "title": "MultilevelEstimators.ML",
    "category": "type",
    "text": "ML()\n\nReturn a multi-level index set.\n\nExamples\n\njulia> ML()\nML\n\nSee also: SL, FT, TD, HC, ZC, AD, U\n\n\n\n\n\n"
},

{
    "location": "manual/#MultilevelEstimators.FT",
    "page": "Manual",
    "title": "MultilevelEstimators.FT",
    "category": "type",
    "text": "FT(d::Integer)\nFT(δ::Real...)\nFT(δ::AbstractVector)\nFT(δ::NTuple)\n\nReturn a full tensor index set in d dimenions with optional weights δ. Default weights are all 1\'s.\n\nExamples\n\njulia> FT(2)\nFT{2}\n\njulia> FT([1, 1/2, 1/2])\nFT{3}\n\njulia> print(FT(2), 4)\n  ◼ ◼ ◼ ◼ ◼\n  ◼ ◼ ◼ ◼ ◼\n  ◼ ◼ ◼ ◼ ◼\n  ◼ ◼ ◼ ◼ ◼\n  ◼ ◼ ◼ ◼ ◼\n\nSee also: SL, ML, TD, HC, ZC, AD, U\n\n\n\n\n\n"
},

{
    "location": "manual/#MultilevelEstimators.TD",
    "page": "Manual",
    "title": "MultilevelEstimators.TD",
    "category": "type",
    "text": "TD(d::Integer)\nTD(δ::Real...)\nTD(δ::AbstractVector)\nTD(δ::NTuple)\n\nExamples\n\njulia> TD(2)\nTD{2}\n\njulia> TD([1, 1/2, 1/2])\nTD{3}\n\njulia> print(TD(2), 4)\n  ◼\n  ◼ ◼\n  ◼ ◼ ◼\n  ◼ ◼ ◼ ◼\n  ◼ ◼ ◼ ◼ ◼\n\nReturn a total degree index set in d dimenions with optional weights δ. Default weights are all 1\'s.\n\nSee also: SL, ML, FT, HC, ZC, AD, U\n\n\n\n\n\n"
},

{
    "location": "manual/#MultilevelEstimators.HC",
    "page": "Manual",
    "title": "MultilevelEstimators.HC",
    "category": "type",
    "text": "HC(d::Integer)\nHC(δ::Real...)\nHC(δ::AbstractVector)\nHC(δ::NTuple)\n\nExamples\n\njulia> HC(2)\nHC{2}\n\njulia> HC([1, 1/2, 1/2])\nHC{3}\n\njulia> print(HC(2), 4)\n  ◼\n  ◼\n  ◼\n  ◼ ◼\n  ◼ ◼ ◼ ◼ ◼\n\nReturn a hyperbolic cross index set in d dimenions with optional weights δ. Default weights are all 1\'s.\n\nSee also: SL, ML, FT, TD, ZC, AD, U\n\n\n\n\n\n"
},

{
    "location": "manual/#MultilevelEstimators.ZC",
    "page": "Manual",
    "title": "MultilevelEstimators.ZC",
    "category": "type",
    "text": "ZC(d::Integer)\nZC(δ::Real...)\nZC(δ::AbstractVector)\nZC(δ::NTuple)\n\nExamples\n\njulia> ZC(2)\nZC{2}\n\njulia> ZC([1, 1/2, 1/2])\nZC{3}\n\njulia> print(ZC(2), 4)\n  ◼ ◼ \n  ◼ ◼ \n  ◼ ◼ ◼ \n  ◼ ◼ ◼ ◼ ◼ \n  ◼ ◼ ◼ ◼ ◼ \n\nReturn a Zaremba cross index set in d dimenions with optional weights δ. Default weights are all 1\'s.\n\nSee also: SL, ML, FT, TD, HC, AD, U\n\n\n\n\n\n"
},

{
    "location": "manual/#MultilevelEstimators.AD",
    "page": "Manual",
    "title": "MultilevelEstimators.AD",
    "category": "type",
    "text": "AD(d::Integer)\n\nReturn an adaptive index set in d dimenions.\n\nExamples\n\njulia> AD(2)\nAD{2}\n\nSee also: SL, ML, FT, TD, HC, ZC, U\n\n\n\n\n\n"
},

{
    "location": "manual/#MultilevelEstimators.U",
    "page": "Manual",
    "title": "MultilevelEstimators.U",
    "category": "type",
    "text": "U(d::Integer)\n\nReturn an unbiased  multi-level or multi-index index set in d dimensions.\n\nExamples\n\njulia> U(2)\nU{2}\n\nSee also: SL, ML, FT, TD, HC, ZC, AD\n\n\n\n\n\n"
},

{
    "location": "manual/#IndexSet-1",
    "page": "Manual",
    "title": "IndexSet",
    "category": "section",
    "text": "SLMLFTTDHCZCADU"
},

]}
