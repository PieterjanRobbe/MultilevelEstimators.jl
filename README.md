# <img src="docs/src/assets/logo.png" alt="alt text" width="75" height="75" align="center"> MultilevelEstimators

| **Documentation** | **Build Status** | **Coverage** |
|-------------------|------------------|--------------|
| [![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://PieterjanRobbe.github.io/MultilevelEstimators.jl/stable) [![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://PieterjanRobbe.github.io/MultilevelEstimators.jl/dev) | [![Build Status](https://travis-ci.org/PieterjanRobbe/MultilevelEstimators.jl.png)](https://travis-ci.org/PieterjanRobbe/MultilevelEstimators.jl) [![Build status](https://ci.appveyor.com/api/projects/status/gh4ka7m9a7qekqu8?svg=true)](https://ci.appveyor.com/project/PieterjanRobbe/multilevelestimators-jl) | [![Coverage](https://codecov.io/gh/PieterjanRobbe/MultilevelEstimators.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PieterjanRobbe/MultilevelEstimators.jl) [![Coverage Status](https://coveralls.io/repos/github/PieterjanRobbe/MultilevelEstimators.jl/badge.svg)](https://coveralls.io/github/PieterjanRobbe/MultilevelEstimators.jl) |

MultilevelEstimators is a Julia package for Multilevel Monte Carlo simulations. It contains the following methods:

+ Monte Carlo (MC) and Quasi-Monte Carlo (QMC)

+ Multilevel Monte Carlo (MLMC) and Multilevel Quasi-Monte Carlo (MLQMC)

+ Multi-Index Monte Carlo (MIMC) and Multi-Index Quasi-Monte Carlo (MIQMC)

+ Adaptive Multi-Index Monte Carlo (A-MIMC) and Adaptive Multi-Index Quasi-Monte Carlo (A-MIQMC)

+ **NEW!** Unbiased Multilevel Monte Carlo (U-MLMC) and Unbiased Multilevel Quasi-Monte Carlo (U-MLQMC)

+ **NEW!** Unbiased Multi-Index Monte Carlo (U-MIMC) and Unbiased Multi-Index Quasi-Monte Carlo (U-MIQMC)


## Installation

The package is not (yet) registered in `METADATA.jl` but can be installed with

```julia
pkg> add https://github.com/PieterjanRobbe/MultilevelEstimators.jl
```

Read the [instructions](https://PieterjanRobbe.github.io/MultilevelEstimators.jl/dev/#Installation-1) in the documentation for more details on how to install MultilevelEstimators and its dependencies.

## Example

An example for a PDE with random coefficients can be found in the [example section](https://PieterjanRobbe.github.io/MultilevelEstimators.jl/dev/example.html#Example-1) in the documentation.

## Documentation

Documentation is available [here](https://PieterjanRobbe.github.io/MultilevelEstimators.jl/dev).

## Related Packages

- [**Reporter.jl**](https://github.com/PieterjanRobbe/Reporter.jl) &mdash; automatic generation of diagnostic information and reports for a MultilevelEstimators simulation
- [**LognormalDiffusionProblems.jl**](https://github.com/PieterjanRobbe/Reporter.jl) &mdash; a package with a complete MultilevelEstimators setup for a 2d elliptic PDE with random coefficients
