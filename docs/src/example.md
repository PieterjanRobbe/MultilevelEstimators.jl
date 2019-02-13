# Example

We consider the elliptic PDE with random coefficients

```math
- \nabla ( a(x, \omega) \nabla u(x, \omega) ) = f(x)
```

defined on the unit square ``D = [0, 1]^2``, with homogeneous boundary conditions ``u(x, \cdot) = 0`` on ``\partial D``. 

Assume the uncertain diffusion coefficient is given as a lognormal random field  ``\log a(x, \omega) = z(x, \omega)`` where ``z(x, \omega)`` is a Gaussian random field (see below).

!!! note

    The following examples assume you have already installed the test dependencies as outlined in the [Installation](@ref) guide.

```@contents
Pages = ["example.md"]
```

## Lognormal diffusion problems

First, we need to define a random field that will be used for the uncertain diffusion coefficient of the PDE. Fortunately, a Julia package that generates and samples from such random fields already exists.

```julia
using GaussianRandomFields
```

To specify the spatial correlation in the random field, we consider the [Matérn](https://en.wikipedia.org/wiki/Matérn_covariance_function) covariance function with a given correlation length and smoothness. This function is predefined in the package. 

```julia
correlation_length = 0.5
smoothness = 2.0
covariance_function = CovarianceFunction(2, Matern(correlation_length, smoothness))
```

We generate samples of this smooth random field using a truncated [Karhunen-Loève expansion](https://en.wikipedia.org/wiki/Karhunen–Loève_theorem) with 250 terms. To this end, we discretize the domain `D` into a regular grid with grid size ``\Delta x = \Delta y = `` 1/128.

```julia
n = 256
pts = range(0, stop=1, length=n)
n_kl = 250
grf = GaussianRandomField(covariance_function, KarhunenLoeve(n_kl), pts, pts)
```

A sample of this Gaussian random field is returned by either calling the `sample`-function directly, or by supplying the random numbers used for sampling from the random field.

```julia
sample(grf)
sample(grf, xi = randn(n_kl))
```

The keyword `xi` corresponds to the sample ``\omega`` in the defintion of the PDE.

We vizualize a sample of this random field.

```julia
using Plots
contourf(sample(grf))
```

```@raw html
<img src="assets/grf.png" width="400">
```

Many more options and details regarding the generation of Gaussian random fields can be found in the [tutorial](https://github.com/PieterjanRobbe/GaussianRandomFields.jl/blob/master/tutorial/tutorial.ipynb) that accomagnies the package [`GaussianRandomFields.jl`](https://github.com/PieterjanRobbe/GaussianRandomFields.jl). 

Suppose we are interested in computing the expected value of a quantity of interest derived from the solution of the PDE. For example, we might be interested in the value of the solution ``u(x, ω)`` at the point ``x=y=``0.5.

## Multilevel Monte Carlo

The basics of any [multilevel method](https://en.wikipedia.org/wiki/Multilevel_Monte_Carlo_method) is a hierarchy of approximations of the model problem with an increasing accuracy, but corresponding increasing cost. For the lognormal diffusion example, such a hierarchy is provided by solving the PDE on an ever finer grid. We define a hierarchy of Gaussian random field generators using
```julia
grfs = Vector{typeof(grf)}(undef, 7)
for i in 1:7
    n = 2^(i+1)
    pts = 1/n:1/n:1-1/n
    grfs[i] = GaussianRandomField(covariance_function, KarhunenLoeve(n_kl), pts, pts)
end
```

The corresponding system matrix that results from a finite difference discretization of the PDE with varying diffusion coefficient is returned by the `elliptic2d`-function from the [`SimpleMultigrid`](https://github.com/PieterjanRobbe/SimpleMultigrid.jl) package.
```julia
using SimpleMultigrid
z = sample(grf)
a = exp.(z)
A = elliptic2d(a)
b = ones(size(A, 1))
u = A\b
contourf(reshape(u, n, n))
```

### Creating an `Estimator`

```julia
using MultilevelEstimators
```

The main type provided by MultilevelEstimators is an [`Estimator`](@ref). This type has two parametric subtypes, `Estimator{<:AbstractIndexSet, <:AbstractSampleMethod}`, where `AbstractIndexSet` is an abstract type for the  index set, and `AbstractSampleMethod` is an abstract type for the sample method. For example, a Multilevel Monte Carlo method will have type signature `Estimator{ML, MC}`.

!!! info
    For a complete list of possible index set types, see [`IndexSet`](@ref IndexSet). For a complete list of possible sample method types, see [`SampleMethod`](@ref SampleMethod).

The main user interaction required by MultilevelEstimators is the definition of a sample function. This function must return a sample of the quantity of interest, and of its difference, for the given discretization parameter (or level) and random parameters. Here is an example for the lognormal diffusion problem:

```julia
function sample_lognormal(level::Level, ω::Vector{<:Real}, grf::GaussianRandomField)

    # solve on finest grid
    z  = sample(grf, xi = ω)
    af = exp.(z)
    Af = elliptic2d(af)
    bf = fill(one(eltype(Af)), size(Af, 1))
    uf = Af\bf
    Qf = uf[length(uf) ÷ 2]

    # compute difference when not on coarsest grid
    dQ = Qf
    if level != Level(0)
        ac = view(af, 2:2:size(af, 1), 2:2:size(af, 2))
        Ac = elliptic2d(ac)
        bc = fill(one(eltype(Af)), size(Ac, 1))
        uc = Ac\bc
        Qc = uc[length(uc) ÷ 2]
        dQ -= Qc
    end
    dQ, Qf
end
```
Thus, when the level parameter is zero, we return the value of the quantity of interest on the coarsest mesh. For a larger level parameter, we return both the difference of the quantity of interest with a coarser grid solution, and the value of the quantity of interest itself.

The sample function can only have two input parameters: a level and a vector of random parameters. For sample functions that require data (such as the precomputed Gaussian random fields in the example above), one might need to add a convenience function:
```julia
sample_lognormal(level, ω) = sample_lognormal(level, ω, grfs[level + one(level)])
```

Finally, for the construction of an `Estimator`, we also need to specify the number of random parameters and their respective distributions. For the lognormal diffusion problem, these random parameters are normally distributed.
```julia
distributions = [Normal() for i in 1:n_kl]
```

!!! info
    A complete list of predefined distributions can be found in the manual for [`Distribution`s](@ref Distribution).

An estimator for the lognormal diffusion problem can thus be created using
```julia
estimator = Estimator(ML(), MC(), sample_lognormal, distributions)
```

### Specifying options

Different options and settings can be passed on to the `Estimator` by supplying the appropriate keyword argument. For example, to set the number of warm up samples to 10, call
```julia
estimator = Estimator(ML(), MC(), sample_lognormal, distributions, nb_of_warm_up_samples = 10)
```

!!! info
    A complete list of optional arguments can be found in the manual for [`Estimator`](@ref).

!!! note

    It is easy to extend the current setup to the computation of the expected value of multiple quantities of interest, by using the keyword `nb_of_qoi`.

### Running a simulation

To start a simulation for the expected value of the quantity of interest up to an absolute (root mean square) error of `5e-3`, call
```julia
h = run(estimator, 5e-3)
```

This function will return a [`History`](@ref) object that contains usefull diagnostics about the simulation.

!!! info
    See [`History`](@ref) for a complete list of diagnostic entries.

!!! note

    By default, the simulation is performed for a sequence of larger tolerances to get improved estimates for the rates of decay of expected value and variance, bias estimation... This sequence is of the form
    ```math
    \text{tol}_i = p^{(n - 1 - i)} \text{tol}, i = 0, \ldots, n - 1
    ``` 
    with ``p>``1 and where the values for ``p`` and ``n`` can be controlled by the keyword arguments `continuation_mul_factor` and `nb_of_tols` respectively. You can disable continuation with the optional argument `continuation = false`. Our experience is that continuation is a very powerful tool when combined with a non-trivial mean square error splitting (`do_mse_splitting = true`) and variance regression (`do_regression = true`). 

!!! note

    By default, samples will be taken in parallel on all available processors (see `addprocs`). The number of processors is controlled by the optional keyword `nb_of_workers`. This can either be a fixed value, or a function, specifying the number of workers to be used on each level.

### Vizualization of the result using [`Reporter`](https://github.com/PieterjanRobbe/Reporter.jl)

MultilevelEstimators can automatically build a set of diagnostic figures based on a [`History`](@ref) file. Load the package by
```julia
using Reporter
```
and generate a report by calling
```julia
report(h)
```

If all goes well, you should now see a webpage with diagnostic information about the simulation, similar to [this one](assets/UntitledEstimator/index.html). Print-ready `.tex`-figures are stored locally under `figures/`.

## Multilevel Quasi-Monte Carlo

[Quasi-Monte Carlo](https://en.wikipedia.org/wiki/Quasi-Monte_Carlo_method) is an alternative way to pick the samples `omega`. The random samples in the Monte Carlo method are replaced by deterministically well-chosen points that increase the accuracy of the estimation with respect to the number of samples. A popular type of such point sets are [rank-1 lattice rules](https://people.cs.kuleuven.be/~dirk.nuyens/qmc-generators/), implemented here as [`LatticeRule32`](@ref). These points are used by default when calling the `QMC`-version of the `Estimator`
```julia
estimator = Estimator(ML(), QMC(), sample_lognormal, distributions)
```
Under the hood, the default lattice rule uses a 3600-dimensional generating vector from [Kuo et al](https://web.maths.unsw.edu.au/~fkuo/lattice/index.html). It is easy to provide your own generating vector with the optional `point_generator`-argument.
```julia
estimator = Estimator(ML(), QMC(), sample_lognormal, point_generator = LatticeRule32("my_gen_vec.txt"))
```

!!! note

    The number of shifts is controlled by the `nb_of_shifts`-keyword. We recommend a value between 10 (default) and 30. The number of samples is updated carefully, with a default multiplication factor of 1.2 (instead of the standard and theoretically justifiable value 2, which we found to be too aggressive). This number can be controlled with the keyword `sample_mul_factor`. 
 
## Multi-Index Monte Carlo

It is easy to extend the lognormal diffusion example to compute multi-index differences. Instead of a difference between a fine and a coarse approximation, we now compute mixed differences.

```julia
function sample_lognormal(index::Index, ω::Vector{<:Real}, grf::GaussianRandomField)

    # solve on finest grid
    z  = sample(grf, xi = ω)
    af = exp.(z)
    Af = elliptic2d(af)
    bf = fill(eltype(Af), size(Af, 1))
    uf = Af\bf
    Qf = uf[length(uf) ÷ 2]

    # compute multi-index differences
    dQ = Qf
    for (key, val) in diff(index)
        step = (index - key).I .+ 1
        ac = view(af, step[1]:step[1]:size(af, 1), step[2]:step[2]:size(af, 2))
        Ac = elliptic2d(ac)
        bc = fill(eltype(Af), size(Ac, 1))
        uc = Ac\bc
        Qc = uc[length(uc) ÷ 2]
        dQ += val * Qc
    end
    dQ, Qf
end
```
All functionality is hidden in the call to `diff(index)`. This function returns a `Dict` with as keys the indices where coarse approximations must be computed, and as values `+1` or `-1`, specifying how these corrections must be added to the quantity of interest on the fine grid.

### Multi-index sets

Different choices for the multi-index sets are available: full tensor index sets ([`FT`](@ref)), total degree index sets ([`TD`](@ref)), hyperbolic cross index sets ([`HC`](@ref)) and Zaremba cross index sets ([`ZC`](@ref)). For example, total degree (TD) index sets look like this 

```@example
using MultilevelEstimators # hide
index_set = TD(2)
for sz = 0:5
    println("size parameter: ", sz)
    print(index_set, sz)
end
```

Weighted index sets are created by specifying appropriate weiths. All weights must be smaller than or equal to 1.
```@example
using MultilevelEstimators # hide
index_set = ZC(1/2, 1)
for sz = 0:5
    println("size parameter: ", sz)
    print(index_set, sz)
end
```

We can create an estimator the usual way:
```julia
estimator = Estimator(TD(2), MC(), sample_lognormal, distributions)
```
where `TD(2)` denotes a 2-dimensional index set.

Before running a simulation, we must precomputed the random fields on all indices.
```julia
grfs = Matrix{typeof(grf)}(undef, 7, 7)
for index in get_index_set(TD(2), 6)
    n = 2 .^(index.I .+ 2)
    pts = map(i -> range(1/i, step = 1/i, stop = 1 - 1/i), n)
    grfs[index + one(index)] = GaussianRandomField(covariance_function, KarhunenLoeve(n_kl), pts...)
end
```

A multi-index Monte Carlo simulation with this estimator is performed by
```julia
h = run(estimator, 5e-3)
```

### Adaptive Multi-index Monte Carlo

Instead of using predefined (weighted or unweighted) index sets, we can also construct the index set adaptively based on some profit indicator (see [6])
```math
P_\ell = \frac{E_\ell}{\sqrt{V_\ell W_\ell}}
```
where ``E_\ell`` and ``V_\ell`` are the expected value and variance of the multi-index difference, respectively, and ``W_\ell`` is the computational cost required to computed a sample of the difference.
```julia
estimator = Estimator(AD(2), MC(), sample_lognormal, distributions)
```

!!! note

    The maximum allowed indices that can be simulated can be set by the keyword `max_search_space`. This keyword should be an index set of one of the above types. The adaptive method will only find indices inside this index set, with the size parameter equal to `max_index_set_param`.
 
!!! note

    Sometimes it is beneficial to add a penalization parameter to the profit indicator to avoid searching too much around the coordinate axes, and bumping into a memory constraint:  
    ```math
    P_\ell = \frac{E_\ell}{(\sqrt{V_\ell W_\ell})^p}
    ```
    with 0 < ``p`` < 1. Specify ``p`` with the optional key `penalization`.
	Another useful optional key is `acceptance_rate`. When this parameter is smaller than 1, a suboptimal index from the active set will be used for further refinement. This index is chosen according to an accept-reject method: we sample a uniform random number using `rand()`, and when this number is larger than the accept ratio, we pich another index at random. Lowering the `acceptance_rate` results in an algorithm with a stronger global search behavior, and is useful when the adaptive method is stuck along one or more directions. Default value is 1 (no globalization). 
 
### Multi-Index Quasi Monte Carlo

Just as in standard Multilevel Monte Carlo, the Multi-Index Monte Carlo method can be extended to use quasi-random numbers instead, see [5]. The setup remains the same, just change `MC()` into `QMC()`. For example, a Hyperbolic Cross Multi-Index Quasi-Monte Carlo estimator can be created by 
```julia
estimator = Estimator(HC(2), QMC(), sample_lognormal, distributions)
```

## Unbiased estimation

The most recent development to MultilevelEstimators is the addition of unbaised multilevel estimators. In contrast to standard Multilevel Monte Carlo, a weighted sum of multilevel contributions is used, thus eliminating the bias in the estimator.

This requires a small change to the sample function, in the sense that it should now return the solutions for all levels samller than or equal to `level`. In the PDE example this can easily be accomplished by using Full Multigrid, see [7] and [8]. See [SimpleMultigrid](https://github.com/PieterjanRobbe/SimpleMultigrid.jl) and [NotSoSimpleMultigrid](https://github.com/PieterjanRobbe/NotSoSimpleMultigrid.jl) for basic example solvers. See also [LognormalDiffusionProblems](https://github.com/PieterjanRobbe/LognormalDiffusionProblems.jl) for a setup using Full Multigrid.

```julia
estimator = estimator(U(1), MC(), sample_lognormal_fmg, distributions)
``` 

!!! warning

    Since this is a recent addition to MultilevelEstimators, it has not been tested as extensively as the other methods. Should you encounter any issues, feel free to submit an issue.

!!! info

    See the index set type [`U`](@ref) for more information and set up.
