# MultilevelEstimators.jl Documentation

For installation instructions, see [Installation](@ref).

For an example on how to use this package, see [Example](@ref).

A full description of the functionality is described in the [Manual](@ref).

## Installation

MultilevelEstimators is not added to the Julia package manager (just yet), but, the package can easily be installed by cloning the git repository. From the Julia REPL, type `]` to enter Pkg mode and run

```julia
pkg> add https://github.com/PieterjanRobbe/MultilevelEstimators.jl.git
```

This will install the main functionality.

!!! note

    The following packages are optional.

For automatic generation of reports and figures you will need the [`Reporter.jl`](https://github.com/PieterjanRobbe/Reporter.jl)  package

```julia
pkg> add https://github.com/PieterjanRobbe/Reporter.jl.git
```

Finally, to run the example problems, you can install the Multigrid solvers and the package to solve lognormal diffusion problems:

```julia
pkg> add https://github.com/PieterjanRobbe/SimpleMultigrid.jl
[...]

pkg> add https://github.com/PieterjanRobbe/NotSoSimpleMultigrid.jl
[...]

pkg> add https://github.com/PieterjanRobbe/LognormalDiffusionProblems.jl
[...]
``` 

## Features

This package features

* Implementation of Multilevel and Multi-Index (Quasi-)Monte Carlo methods.

* Native support for parallelization on multicore machines.

* Full control of advanced algorithm options such as variance regression, mean square error splitting, ...

* Automatic generation of reports with print-quality figures using [`Reporter.jl`](https://github.com/PieterjanRobbe/Reporter.jl).

Most recent addition to the package is **unbiased** Multilevel and Multi-Index Monte Carlo with automatic learning of the discrete distribution of samples across all levels or mult-indices. (See index set type [`U`](@ref).)

## References

The algorithms implemented in this package are loosely based on the following papers:

Basic Multilevel Monte Carlo method:

* Giles, M. B. *Multilevel Monte Carlo Path Simulation*. Operations Research 56.3 (2008): 607-617.

Variance regression and continuation:

* Collier, N., Haji-Ali, A. L., Nobile, F., von Schwerin, E., and Tempone, R. *A Continuation Multilevel Monte Carlo Algorithm*. BIT Numerical Mathematics 55.2 (2015): 399-432. 

Quasi-Monte Carlo and Multilevel Quasi-Monte Carlo methods:

* Giles, M. B., and Waterhouse, B. J. *Multilevel Quasi-Monte Carlo Path Simulation*. Advanced Financial Modelling, Radon Series on Computational and Applied Mathematics (2009): 165-181.

Multi-Index (Quasi-)Monte Carlo methods:

* Robbe, P., Nuyens, D., and Vandewalle, S. *A Multi-Index Quasi-Monte Carlo Algorithm for Lognormal Diffusion Problems*. SIAM Journal on Scientific Computing 39.5 (2017): S851-S872.

Adaptive Multi-Index methods:

* Robbe, P., Nuyens, D., and Vandewalle, S. *A Dimension-Adaptive Multi-Index Monte Carlo Method Applied to a Model of a Heat Exchanger*. International Conference on Monte Carlo and Quasi-Monte Carlo Methods in Scientific Computing. Springer Proceedings in Mathematics & Statistics 241 (2018): 429-445.

Unbiased estimation:

* Robbe, P., Nuyens, D. and Vandewalle, S. *Recycling Samples in the Multigrid Multilevel (Quasi-)Monte Carlo Method*. SIAM Journal on Scientific Computing, to appear (2019).

* Robbe, P., Nuyens, D. and Vandewalle, S. *Enhanced Multi-Index Monte Carlo by means of Multiple Semi-coarsened Multigrid for Anisotropic Diffusion Problems*. In preparation (2019).
