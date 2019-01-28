# Example

As a benchmark problem, we consider the elliptic PDE with random coefficients

```math
- \nabla ( a(x, \omega) \nabla u(x, \omega) ) = f(x)
```

defined on the unit square ``D = [0, 1]^2``, with homogeneous boundary conditions ``u(x, \cdot) = 0`` on ``\partial D``. 

The uncertain diffusion coefficient is given as a lognormal random field  ``\log a(x, \omega) = z(x, \omega)`` where ``z(x, \omega)`` is a Gaussian random field (see below).

!!! note

    The following examples assume you have already installed the test dependencies as outlined in the [Installation](@ref) guide.


## Lognormal diffusion problems

First, we need to define a random field that will be used for the uncertain diffusion coefficient of the PDE. Fortunately, a Julia package that generates and samples from such random fields already exists.

```julia
using GaussianRandomFields
```

To specify the spatial correlation in the random field, we consider the [Mat\'ern](https://en.wikipedia.org/wiki/Matérn_covariance_function) covariance function with a given correlation length and smoothness. This function is predefined in the package. 

```julia
corr_length = 0.5
smoothness = 2.0
cov_fun = CovarianceFunction(2, Matern(corr_length, smoothness))
```

We generate samples of this smooth random field using a [Karhunen-Lo\`eve expansion](https://en.wikipedia.org/wiki/Karhunen–Loève_theorem) with 250 terms. To this end, we discretize the domain `D` into a regular grid with grid size ``\Delta x = \Delta y = `` 1/128.

```julia
n = 128
pts = range(0, stop=1, length=n)
grf = GaussianRandomField(cov_fun, KarhunenLoeve(250), pts, pts)
```

A sample of this Gaussian random field is returned by either calling the `sample`-function directly, or by supplying the random numbers used for sampling from the random field.

```julia
sample(grf)
sample(grf, xi = randn(250))
```

The keyword `xi` corresponds to the sample ``\omega`` in the defintion of the PDE.

We vizualize a sample of this random field.

```julia
using Plots
contourf(sample(grf))
```

```@raw html
<img src="../assets/grf.png" width="400">
```

Many more options and details regarding the generation of Gaussian random fields can be found in the [tutorial](https://github.com/PieterjanRobbe/GaussianRandomFields.jl/blob/master/tutorial/tutorial.ipynb) that accomagnies the package [`GaussianRandomFields.jl`](https://github.com/PieterjanRobbe/GaussianRandomFields.jl). 

Suppose we are interested in computing the expected value of a quantity of interest derived from the solution of the PDE. For example, we might be interested in the value of the solution ``u(x, \omega)`` at the point ``x=y=0.5``.

## Multilevel Monte Carlo

The basics of any [multilevel method](https://en.wikipedia.org/wiki/Multilevel_monte_carlo_method) is a hierarchy of approximations of the model problem with an increasing accuracy, but corresponding increasing cost. For the lognormal diffusion example, such a hierarchy is provided by solving the PDE on an ever finer grid. We define a hierarchy of Gaussian random field generators using
```julia
grfs = [GaussianRandomField(C, KarhunenLoeve(250), 1/i:1/i:1-1/i, 1/i:1/i:1-1/i) for i in 2 .^(2:8)]
```

The corresponding system matrix that results from a finite difference discretization of the PDE is returned by the `elliptic2d`-function from the [`SimpleMultigrid`](https://github.com/PieterjanRobbe/SimpleMultigrid.jl) package.
```julia
z = sample(grf)
a = exp(z)
A = elliptic2d(a)
```

### Creating an `Estimator`

The main type provided by `MultilevelEstimators.jl` is an [`Estimator`](@ref). This type has two parametric subtypes, `Estimator{<:AbstractIndexSet, <:AbstractSampleMethod}`, where `AbstractIndexSet` is an abstract type for the  index set, and `AbstractSampleMethod` is an abstract type for the sample method. For example, a Multilevel Monte Carlo method will have type signature `Estimator{ML, MC}`. For a complete list of possible index set types, see [`IndexSet`](@ref). For a complete list of possible sample method types, see [`SampleMethod`](@ref).

The main user interaction required by `MultilevelEstimators.jl` is the definition of a sample function. This function must return a sample of the quantity of interest, and of its difference, for the given discretization parameter (or level) and random parameters. Here is an example for the lognormal diffusion problem:

```julia
function sample_lognormal(level::Level, \omega::Vector{<:Real}, grf::GaussianRandomField)

    # solve on finest grid
    z  = sample(grf, xi=\omega)
    af = exp(z)
    Af = elliptic2d(af)
    bf = fill(eltype(Af), size(Af, 1))
    uf = Af\bf
    Qf = uf[length(uf) >> 1]

    if level == 0
        Qc = zero(Qf)
    else
        ac = view(af, 2:2:end, 2:2:end)
        Ac = elliptic2d(ac)
        bc = fill(eltype(Af), size(Ac, 1))
        uc = Ac\bc
        Qc = uc[length(uc) >> 1]
    end
    Qf - Qc
end
```
When the level parameter is zero, we return the value of the quantity of interest on the coarsest mesh. For a larger level parameter, we return both the difference of the quantity of interest with a coarser grid solution, and the value of the quantity of interest itself.

The sample function can only have two input parameters: a level and a vector of random parameters. For sample functions that require data (such as the precomputed Gaussian random fields in the example above), one might need to add a convenience function:
```julia
sample_lognormal(level, \omega) = sample_lognormal(level, \omega, grfs[level + one(level)])
```

Finally, for the construction of an `Estimator`, we also need to specify the number of random parameters and their distribution. For the lognormal diffusion problem, these random parameters are normally distributed.
```julia
distributions = [Normal() for i in 1:250]
```
A full list of predefined distributions can be found in the manual for [`Distributions`](@ref).

An estimator for the lognormal diffusion problem can thus be created using
```julia
estimator = Estimator(ML(), MC(), sample_lognormal, distributions)
```

### Specifying options

Different options and settings can be passed on to the `Estimator` by supplying the appropriate keyword argument. For example, to set the number of warm up samples to 10, call
```julia
estimator = Estimator(ML(), MC(), sample_lognormal, distributions, nb_of_warm_up_samples = 10)
```
A complete list of optional arguments can be found in the manual for [`Estimator`](@ref).

### Running a simulation

To start a simulation that computes the expected value of the quantity of interest up to an absolute (root mean square) error of `1e-3`, call
```julia
h = run(estimator, 1e-3)
```

This function will return a [`History`](@ref) object that contains usefull diagnostics about the simulation.
See [`History`](@ref) for a complete list of diagnostic entries.

### Vizualization of the result using [`Reporter`](@ref)

`MultilevelEstimators.jl` can automatically build a set of diagnostic figures based on a [`History`](@ref) file. Load the package by
```julia
using Reporter
```
and generate a report by calling
```julia
report(h)
```

If all goes well, you should now see a webpage with diagnostic information about the simulation, similar to [this one](somelink). Print-ready `.tex`-figures are stored locally under `figures/`.

## Multilevel Quasi-Monte Carlo

[Quasi-Monte Carlo](wiki) is an alternative way to pick the samples `omega`. The random samples in the Monte Carlo method are replaced by deterministically well-chosen points that increase the accuracy of the estimation with respect to the number of samples. A well-known type of such point sets are [`rank-1 lattice rules](), implemented here as [`LatticeRule32`](@ref). These points are used by default when calling the `QMC`-version of the `Estimator`
```julia
estimator = Estimator(ML(), QMC(), sample_lognormal, distributions)
```
Under the hood, the default lattice rule uses a 3600-dimensional generating vector from [Kuo et al.](). It is easy to provide your own generating vector with the optional `point_generator`-argument.
```julia
estimator = Estimator(ML(), QMC(), sample_lognormal, point_generator=LatticeRule32("my_gen_vec.txt"))
```

The number of shifts is controlled by the `nb_of_shifts`-keyword. We recommand a value between 10 (default) and 30. The number of samples is updated carefully, with a default multiplication factor of 1.2 (instead of the standard and theoretically justifiable value 2, which we found to be too aggressive). This number can be controlled with the keyword `sample_mul_factor`. 
 
## Multi-Index Monte Carlo

### Multi-index sets
```@example
using MultilevelEstimators # hide
index_set = TD(2)
for L = 0:5
    println("level: ", L)
    print(index_set, L)
end
```

lognormal example

### Adaptive Multi-index Monte Carlo
with penalization

## Unbiased estimation

## Parallelization

