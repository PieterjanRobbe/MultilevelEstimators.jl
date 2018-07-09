# Multilevel Monte Carlo for Julia
[![Build Status](https://travis-ci.org/PieterjanRobbe/MultilevelEstimators.jl.png)](https://travis-ci.org/PieterjanRobbe/MultilevelEstimators.jl)

This module provides a versatile implementation of Multilevel Monte Carlo (MLMC) methods and Multi-Index Monte Carlo (MIMC) methods, and their Quasi-Monte Carlo (QMC) counterparts, Multilevel Quasi-Monte Carlo (MLQMC) and Multi-Index Quasi-Monte Carlo (MIQMC). The module is mainly aimed at solving partial differential equations (PDEs) with random coefficients, but can be used for any forward uncertainty quantification (UQ) problem by specifying the appropriate sample function.

The module has been tested for challenging 3D problems with up to 5000 uncertainties, and the parallel implementation was tested on a compute node with 24 cores.

Recent additions are Adaptive Multi-Index Monte Carlo (AMIMC) and Multigrid Multilevel Monte Carlo (MG-MLMC).

## Installation

```julia
Pkg.clone("https://github.com/PieterjanRobbe/MultilevelEstimators.jl")
```

Load the package in Julia by

```julia
using MultilevelEstimators
```

## Usage

### Introduction to Multilevel Monte Carlo

In most mathematical models, parameters or coefficients are unknown or subject to uncertainty, particularly due to lack of data or measurements. Often, these problems involve the computation of a *quantity of interest* as the expected value over the uncertain input parameters. The classical sample-based approach then chooses <img src="https://latex.codecogs.com/svg.latex?\Large&space;N"/> realisations of the uncertain parameters and approximates this expected value as a sample average. The Multilevel Monte Carlo (MLMC) method improves the error versus work complexity rate of the classical approach by using models  with different levels of accuracy. These models are called *levels*. The main idea is to write the approximation to <img src="https://latex.codecogs.com/svg.latex?\Large&space;G"/> at the most accurate level <img src="https://latex.codecogs.com/svg.latex?\Large&space;L"/> as a telescoping sum

<center>
<img src="https://latex.codecogs.com/svg.latex?\Large&space;\mathbb{E}[G_L(\omega)]=\mathbb{E}[G_0(\omega)]+\sum_{\ell=1}^L\mathbb{E}[G_\ell(\omega)-G_{\ell-1}(\omega)]=\sum_{\ell=0}^L\mathbb{E}[\Delta\,G_\ell(\omega)]"/>
</center>

Note that <img src="https://latex.codecogs.com/svg.latex?\Large&space;\ell=0"/> is the least accurate model. Hence, instead of approximating the expected value of the quantity of interest on the finest mesh, the MLMC method approximates differences <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Delta\,G_\ell"/> at different levels <img src="https://latex.codecogs.com/svg.latex?\Large&space;\ell"/>. If the variance of these differences goes sufficiently fast to zero as <img src="https://latex.codecogs.com/svg.latex?\Large&space;\ell"/> increases, most samples are taken at models with low accuracy, hence low cost. Typically, only very few samples are needed at the finest mesh.


### A First Simple MLMC Example

PDEs with random coefficients typically describe cell movements, fluid pressures, temperatures or flow rates, all subject to uncertainty. For example, a model for fluid flow through porous media, commonly used in geophysics, is the elliptic PDE

<img src="https://latex.codecogs.com/svg.latex?\Large&space;-\nabla\cdot(k(x,\omega)\nabla\,p(x,\omega))=f(x)"/>

with random diffusion coefficient <img src="https://latex.codecogs.com/svg.latex?\Large&space;k(x,\omega)"/>. Here, <img src="https://latex.codecogs.com/svg.latex?\Large&space;k"/> represents the uncertain permeability of the porous material, and <img src="https://latex.codecogs.com/svg.latex?\Large&space;p(x,\omega)"/> is the unknown pressure head. The problem is physically meaningful if <img src="https://latex.codecogs.com/svg.latex?\Large&space;k(x,\omega)"/>> 0, hence, usually <img src="https://latex.codecogs.com/svg.latex?\Large&space;k(x,\omega)=\exp(Z(x,\omega))"/>, where <img src="https://latex.codecogs.com/svg.latex?\Large&space;Z"/> is an underlying Gaussian random field with known characteristics. This is known as the *lognormal* case.

We generate samples of the Gaussian random field <img src="https://latex.codecogs.com/svg.latex?\Large&space;Z(x,\omega)"/> using the Julia-package [GaussianRandomFields.jl](https://github.com/PieterjanRobbe/GaussianRandomFields.jl). Below are some typical examples:

<p align="center">
<img src="https://github.com/PieterjanRobbe/GaussianRandomFields.jl/blob/master/figures/front.png" width="650" align="middle">
</p>

When dealing with PDEs, the natural sequence of models with different accuracy, required in the Multilevel Monte Carlo method, can be obtained by varying the discretization mesh. For example, one can use a Finite Difference (FD) method and consider meshes with a decreasing mesh size parameter as the level `l` increases. Let us impose a geometrical grid hierarchy as follows:

<center>

| Level         | Number of DOF's           |
| ------------- |:-------------:|
| `0`           | `4`x`4` |
| `1`           | `8`x`8` |
| `2`           | `16`x`16` |
|   ...          | ... |
| `l`          | `(4*2^l)`x`(4*2^l)` |

</center>

This example is already defined in the [applications](https://github.com/PieterjanRobbe/MultilevelEstimators.jl/tree/fresh/applications/SPDE) folder. More details and an example to study can be found there. You can load the example in Julia by entering

```julia
push!(LOAD_PATH,Pkg.dir(joinpath("MultilevelEstimators","applications","SPDE")))
using SPDE
```

Initialize the MLMC estimator by calling the appropriate `init` function

```julia
srand(2018) # for reproducability
estimator = init_lognormal_diffusion_mlmc()
```

The quantity of interest in this case is a point evaluation of the pressure in the middle of the domain. We can compute the expected value of this quantity up to a tolerance of `tol = 1e-3` by

```julia
h = run(estimator,1e-3)
```

The output of this simulation is

```julia
--------------------------------------------------------------------------------
*** MultilevelEstimators.jl @2018-07-09T15:19:33.006
*** This is a Multilevel Monte Carlo estimator
*** Simulating SPDE.lognormal_diffusion_single
*** Tolerance on RMSE ϵ = 1.000e-03
--------------------------------------------------------------------------------
Taking 20 warm-up samples at level 0...
--------------------------------------------------------------------------------
  level       E              V              N              W              
--------------------------------------------------------------------------------
  (0,)        6.72944e-02    3.25352e-03    20             8.00000e+00    
--------------------------------------------------------------------------------
Samples will be updated according to
-----------------------------
  level       N              
-----------------------------
  (0,)        6508
-----------------------------
Taking 6488 additional samples at level 0...
--------------------------------------------------------------------------------
  level       E              V              N              W              
--------------------------------------------------------------------------------
  (0,)        7.68503e-02    6.09148e-03    6508           8.00000e+00    
--------------------------------------------------------------------------------
Checking convergence...
  ==> Rates: α ≈   NaN, β ≈   NaN, γ ≈   NaN.
  ==> Variance of the estimator ≈ 9.35998e-07.
  ==> Bias of the estimator ≈         NaN.
Taking 20 warm-up samples at level 1...
--------------------------------------------------------------------------------
  level       E              V              N              W              
--------------------------------------------------------------------------------
  (0,)        7.68503e-02    6.09148e-03    6508           8.00000e+00    
  (1,)        7.11509e-03    8.66849e-05    20             3.06274e+01    
--------------------------------------------------------------------------------
Samples will be updated according to
-----------------------------
  level       N              
-----------------------------
  (0,)        15027
  (1,)        917
-----------------------------
Taking 8519 additional samples at level 0...
Taking 897 additional samples at level 1...
--------------------------------------------------------------------------------
  level       E              V              N              W              
--------------------------------------------------------------------------------
  (0,)        7.74528e-02    6.23115e-03    15027          8.00000e+00    
  (1,)        9.34629e-03    2.12659e-04    917            3.06274e+01    
--------------------------------------------------------------------------------
Checking convergence...
  ==> Rates: α ≈   NaN, β ≈   NaN, γ ≈   NaN.
  ==> Variance of the estimator ≈ 6.46571e-07.
  ==> Bias of the estimator ≈         NaN.
Taking 20 warm-up samples at level 2...
--------------------------------------------------------------------------------
  level       E              V              N              W              
--------------------------------------------------------------------------------
  (0,)        7.74528e-02    6.23115e-03    15027          8.00000e+00    
  (1,)        9.34629e-03    2.12659e-04    917            3.06274e+01    
  (2,)        2.58469e-03    1.16025e-05    20             6.41966e+01    
--------------------------------------------------------------------------------
Samples will be updated according to
-----------------------------
  level       N              
-----------------------------
  (0,)        18491
  (1,)        1746
  (2,)        282
-----------------------------
Taking 3464 additional samples at level 0...
Taking 829 additional samples at level 1...
Taking 262 additional samples at level 2...
--------------------------------------------------------------------------------
  level       E              V              N              W              
--------------------------------------------------------------------------------
  (0,)        7.72633e-02    6.04055e-03    18491          8.00000e+00    
  (1,)        8.80863e-03    2.33218e-04    1746           3.06274e+01    
  (2,)        2.75612e-03    1.32898e-05    282            6.41966e+01    
--------------------------------------------------------------------------------
Checking convergence...
  ==> Rates: α ≈ 1.676, β ≈ 4.133, γ ≈ 1.068.
  ==> Variance of the estimator ≈ 5.07375e-07.
  ==> Bias of the estimator ≈ 8.62360e-04.
No convergence yet. RMSE ≈ 1.11850e-03 > 1.000e-03.
Adding an extra level...
Taking 20 warm-up samples at level 3...
--------------------------------------------------------------------------------
  level       E              V              N              W              
--------------------------------------------------------------------------------
  (0,)        7.72633e-02    6.04055e-03    18491          8.00000e+00    
  (1,)        8.80863e-03    2.33218e-04    1746           3.06274e+01    
  (2,)        2.75612e-03    1.32898e-05    282            6.41966e+01    
  (3,)        1.03354e-03    9.46216e-07    20             1.05569e+02    
--------------------------------------------------------------------------------
Samples will be updated according to
-----------------------------
  level       N              
-----------------------------
  (0,)        11109
  (1,)        1116
  (2,)        184
  (3,)        39
-----------------------------
Taking 19 additional samples at level 3...
--------------------------------------------------------------------------------
  level       E              V              N              W              
--------------------------------------------------------------------------------
  (0,)        7.72633e-02    6.04055e-03    18491          8.00000e+00    
  (1,)        8.80863e-03    2.33218e-04    1746           3.06274e+01    
  (2,)        2.75612e-03    1.32898e-05    282            6.41966e+01    
  (3,)        8.92268e-04    6.14640e-07    39             1.05569e+02    
--------------------------------------------------------------------------------
Checking convergence...
  ==> Rates: α ≈ 1.652, β ≈ 4.284, γ ≈ 0.893.
  ==> Variance of the estimator ≈ 5.23135e-07.
  ==> Bias of the estimator ≈ 2.88863e-04.
  ==> Non-trivial MSE splitting parameter ≈ 0.85.
--------------------------------------------------------------------------------
  level       E              V              N              W              
--------------------------------------------------------------------------------
  (0,)        7.72633e-02    6.04055e-03    18491          8.00000e+00    
  (1,)        8.80863e-03    2.33218e-04    1746           3.06274e+01    
  (2,)        2.75612e-03    1.32898e-05    282            6.41966e+01    
  (3,)        8.92268e-04    6.14640e-07    39             1.05569e+02    
--------------------------------------------------------------------------------
Convergence reached. RMSE ≈ 7.78831e-04.
--------------------------------------------------------------------------------
*** MultilevelEstimators.jl @2018-07-09T15:20:06.993
*** Successfull termination
--------------------------------------------------------------------------------
```

The MLMC algorithm used four levels and took 18491 PDE solves on a grid with 16 unknowns, 1746 PDE solves on a grid with 64 unknowns, 282 PDE solves on a grid with 256 unknowns, and 39 PDE solves on a grid with 1024 unkowns. A call to `run` runs the imposed hierarchy and returns a `History`-object containg various log-data. For example, the result of the simulation is 0.0897 with a variance of 0.00628.

```julia
h[:mean]
h[:var]
```

Note that additional calls to `run` will reuse the samples already taken on previous runs.
A full list of options that can be passed to the estimator is given below. For optimal efficiency, set `continuate = true`, `do_regression = true` and `do_splitting = true`. We also provide a `reset` function to clear the samples already taken:

```julia
clear(estimator)
```

### Multi-Index Monte Carlo

Multi-Index Monte Carlo (MIMC) is a generalisation of MLMC where the levels are replaced by multi-indices. Whereas MLMC uses levels that refine in both <img src="https://latex.codecogs.com/svg.latex?\Large&space;x"/>- and <img src="https://latex.codecogs.com/svg.latex?\Large&space;y"/>-direction simultaneously, the MIMC method will allow for grids that refine only in <img src="https://latex.codecogs.com/svg.latex?\Large&space;x"/> or only in <img src="https://latex.codecogs.com/svg.latex?\Large&space;y"/>. For non-isotropic examples, where the MLMC method would perform badly, MIMC allows to again achieve the optimal convergence rate of error versus work. An illustration of these grids is shown below:

<p align="center">
<img align="middle" src="https://github.com/PieterjanRobbe/MultilevelEstimators.jl/blob/master/figures/fig_mimc.png" width="450">
</p>

The classical MLMC method would only consider meshes on the main diagonal, whereas MIMC considers a subset of all grids depicted. This subset is called the *index set*. As for the Smolyak construction in Sparse Grids, it can be shown that an index set consisting of the grids in the left upper triangle (the `TD` case) are, in some sense, optimal, see [2]. We provide 6 different index sets, `SL` (single level), `ML` (multilevel), `FT` (full tensor), `TD` (total degree), `HC` (hyperbolic cross) and `AD`. This last index set is an Adaptive Multi-Index Monte Carlo (AMIMC) method, see below.

Let us solve the model problem from the previous section using MIMC with TD index sets.

```julia
srand(2018) # for reproducability
estimator = init_lognormal_diffusion_mimc()
h = run(estimator,1e-3)
```

The output of the simulation is

```julia
--------------------------------------------------------------------------------
*** MultilevelEstimators.jl @2018-07-09T15:30:50.76
*** This is a Multi-Index Monte Carlo estimator (TD index set)
*** Simulating SPDE.lognormal_diffusion_single
*** Tolerance on RMSE ϵ = 1.000e-03
--------------------------------------------------------------------------------
Taking 20 warm-up samples at index (0, 0)...
Shape of the index set:
  ◼ 
--------------------------------------------------------------------------------
  index       E              V              N              W              
--------------------------------------------------------------------------------
  (0, 0)      6.72944e-02    3.25352e-03    20             8.00000e+00    
--------------------------------------------------------------------------------
Samples will be updated according to
-----------------------------
  index       N              
-----------------------------
  (0, 0)      6508
-----------------------------
Taking 6488 additional samples at index (0, 0)...
Shape of the index set:
  ◼ 
--------------------------------------------------------------------------------
  index       E              V              N              W              
--------------------------------------------------------------------------------
  (0, 0)      7.68503e-02    6.09148e-03    6508           8.00000e+00    
--------------------------------------------------------------------------------
Checking convergence...
  ==> Rates: α ≈ (   NaN,   NaN), β ≈ (   NaN,   NaN), γ ≈ (   NaN,   NaN).
  ==> Variance of the estimator ≈ 9.35998e-07.
  ==> Bias of the estimator ≈         NaN.
Taking 20 warm-up samples at index (1, 0)...
Taking 20 warm-up samples at index (0, 1)...
Shape of the index set:
  ◼ 
  ◼ ◼ 
--------------------------------------------------------------------------------
  index       E              V              N              W              
--------------------------------------------------------------------------------
  (0, 0)      7.68503e-02    6.09148e-03    6508           8.00000e+00    
  (0, 1)      5.82821e-03    5.94196e-05    20             4.66274e+01    
  (1, 0)      4.61627e-03    2.54315e-04    20             4.66274e+01    
--------------------------------------------------------------------------------
Samples will be updated according to
-----------------------------
  index       N              
-----------------------------
  (0, 0)      21098
  (0, 1)      864
  (1, 0)      1786
-----------------------------
Taking 14590 additional samples at index (0, 0)...
Taking 844 additional samples at index (0, 1)...
Taking 1766 additional samples at index (1, 0)...
Shape of the index set:
  ◼ 
  ◼ ◼ 
--------------------------------------------------------------------------------
  index       E              V              N              W              
--------------------------------------------------------------------------------
  (0, 0)      7.72061e-02    5.93900e-03    21098          8.00000e+00    
  (0, 1)      5.23463e-03    1.27874e-04    864            4.66274e+01    
  (1, 0)      4.69202e-03    1.26412e-04    1786           4.66274e+01    
--------------------------------------------------------------------------------
Checking convergence...
  ==> Rates: α ≈ (   NaN,   NaN), β ≈ (   NaN,   NaN), γ ≈ (   NaN,   NaN).
  ==> Variance of the estimator ≈ 5.00277e-07.
  ==> Bias of the estimator ≈         NaN.
Taking 20 warm-up samples at index (2, 0)...
Taking 20 warm-up samples at index (1, 1)...
Taking 20 warm-up samples at index (0, 2)...
Shape of the index set:
  ◼ 
  ◼ ◼ 
  ◼ ◼ ◼ 
--------------------------------------------------------------------------------
  index       E              V              N              W              
--------------------------------------------------------------------------------
  (0, 0)      7.72061e-02    5.93900e-03    21098          8.00000e+00    
  (0, 1)      5.23463e-03    1.27874e-04    864            4.66274e+01    
  (0, 2)      1.45688e-03    7.87958e-06    20             8.01966e+01    
  (1, 0)      4.69202e-03    1.26412e-04    1786           4.66274e+01    
  (1, 1)      1.08691e-03    1.32611e-05    20             9.98823e+01    
  (2, 0)      5.88808e-04    1.51060e-05    20             8.01966e+01    
--------------------------------------------------------------------------------
Samples will be updated according to
-----------------------------
  index       N              
-----------------------------
  (0, 0)      25520
  (0, 1)      1552
  (0, 2)      294
  (1, 0)      1543
  (1, 1)      342
  (2, 0)      407
-----------------------------
Taking 4422 additional samples at index (0, 0)...
Taking 688 additional samples at index (0, 1)...
Taking 274 additional samples at index (0, 2)...
Taking 322 additional samples at index (1, 1)...
Taking 387 additional samples at index (2, 0)...
Shape of the index set:
  ◼ 
  ◼ ◼ 
  ◼ ◼ ◼ 
--------------------------------------------------------------------------------
  index       E              V              N              W              
--------------------------------------------------------------------------------
  (0, 0)      7.71469e-02    6.01058e-03    25520          8.00000e+00    
  (0, 1)      5.45626e-03    1.24835e-04    1552           4.66274e+01    
  (0, 2)      1.59881e-03    9.36213e-06    294            8.01966e+01    
  (1, 0)      4.69202e-03    1.26412e-04    1786           4.66274e+01    
  (1, 1)      4.01388e-04    1.72171e-05    342            9.98823e+01    
  (2, 0)      1.01867e-03    1.43615e-05    407            8.01966e+01    
--------------------------------------------------------------------------------
Checking convergence...
  ==> Rates: α ≈ ( 2.204, 1.771), β ≈ ( 3.138, 3.737), γ ≈ ( 0.782, 0.782).
  ==> Variance of the estimator ≈ 5.04211e-07.
  ==> Bias of the estimator ≈ 8.98042e-04.
No convergence yet. RMSE ≈ 1.14485e-03 > 1.000e-03.
Adding an extra level...
Taking 20 warm-up samples at index (3, 0)...
Taking 20 warm-up samples at index (2, 1)...
Taking 20 warm-up samples at index (1, 2)...
Taking 20 warm-up samples at index (0, 3)...
Shape of the index set:
  ◼ 
  ◼ ◼ 
  ◼ ◼ ◼ 
  ◼ ◼ ◼ ◼ 
--------------------------------------------------------------------------------
  index       E              V              N              W              
--------------------------------------------------------------------------------
  (0, 0)      7.71469e-02    6.01058e-03    25520          8.00000e+00    
  (0, 1)      5.45626e-03    1.24835e-04    1552           4.66274e+01    
  (0, 2)      1.59881e-03    9.36213e-06    294            8.01966e+01    
  (0, 3)      3.17754e-04    1.17917e-06    20             1.21569e+02    
  (1, 0)      4.69202e-03    1.26412e-04    1786           4.66274e+01    
  (1, 1)      4.01388e-04    1.72171e-05    342            9.98823e+01    
  (1, 2)      1.32163e-04    1.62660e-06    20             1.89648e+02    
  (2, 0)      1.01867e-03    1.43615e-05    407            8.01966e+01    
  (2, 1)     -5.59629e-04    5.36879e-06    20             1.89648e+02    
  (3, 0)      5.07976e-04    7.47523e-07    20             1.21569e+02    
--------------------------------------------------------------------------------
Samples will be updated according to
-----------------------------
  index       N              
-----------------------------
  (0, 0)      15121
  (0, 1)      903
  (0, 2)      189
  (0, 3)      55
  (1, 0)      909
  (1, 1)      230
  (1, 2)      52
  (2, 0)      234
  (2, 1)      93
  (3, 0)      44
-----------------------------
Taking 35 additional samples at index (0, 3)...
Taking 32 additional samples at index (1, 2)...
Taking 73 additional samples at index (2, 1)...
Taking 24 additional samples at index (3, 0)...
Shape of the index set:
  ◼ 
  ◼ ◼ 
  ◼ ◼ ◼ 
  ◼ ◼ ◼ ◼ 
--------------------------------------------------------------------------------
  index       E              V              N              W              
--------------------------------------------------------------------------------
  (0, 0)      7.71469e-02    6.01058e-03    25520          8.00000e+00    
  (0, 1)      5.45626e-03    1.24835e-04    1552           4.66274e+01    
  (0, 2)      1.59881e-03    9.36213e-06    294            8.01966e+01    
  (0, 3)      3.94629e-04    1.37699e-06    55             1.21569e+02    
  (1, 0)      4.69202e-03    1.26412e-04    1786           4.66274e+01    
  (1, 1)      4.01388e-04    1.72171e-05    342            9.98823e+01    
  (1, 2)      3.71712e-04    3.72400e-06    52             1.89648e+02    
  (2, 0)      1.01867e-03    1.43615e-05    407            8.01966e+01    
  (2, 1)      1.52900e-04    4.70600e-06    93             1.89648e+02    
  (3, 0)      4.38346e-04    6.49751e-07    44             1.21569e+02    
--------------------------------------------------------------------------------
Checking convergence...
  ==> Rates: α ≈ ( 1.710, 1.895), β ≈ ( 3.802, 3.251), γ ≈ ( 0.691, 0.691).
  ==> Variance of the estimator ≈ 6.66232e-07.
  ==> Bias of the estimator ≈ 6.10507e-04.
  ==> Non-trivial MSE splitting parameter ≈ 0.99.
No convergence yet. RMSE ≈ 1.01929e-03 > 1.000e-03.
Adding an extra level...
Taking 11 warm-up samples at index (4, 0)...
Taking 20 warm-up samples at index (3, 1)...
Taking 15 warm-up samples at index (2, 2)...
Taking 20 warm-up samples at index (1, 3)...
Taking 16 warm-up samples at index (0, 4)...
Shape of the index set:
  ◼ 
  ◼ ◼ 
  ◼ ◼ ◼ 
  ◼ ◼ ◼ ◼ 
  ◼ ◼ ◼ ◼ ◼ 
--------------------------------------------------------------------------------
  index       E              V              N              W              
--------------------------------------------------------------------------------
  (0, 0)      7.71469e-02    6.01058e-03    25520          8.00000e+00    
  (0, 1)      5.45626e-03    1.24835e-04    1552           4.66274e+01    
  (0, 2)      1.59881e-03    9.36213e-06    294            8.01966e+01    
  (0, 3)      3.94629e-04    1.37699e-06    55             1.21569e+02    
  (0, 4)      9.18491e-05    1.89814e-08    16             1.69443e+02    
  (1, 0)      4.69202e-03    1.26412e-04    1786           4.66274e+01    
  (1, 1)      4.01388e-04    1.72171e-05    342            9.98823e+01    
  (1, 2)      3.71712e-04    3.72400e-06    52             1.89648e+02    
  (1, 3)      5.21324e-05    3.23617e-07    20             2.72393e+02    
  (2, 0)      1.01867e-03    1.43615e-05    407            8.01966e+01    
  (2, 1)      1.52900e-04    4.70600e-06    93             1.89648e+02    
  (2, 2)     -2.74190e-04    1.13364e-06    15             2.15217e+02    
  (3, 0)      4.38346e-04    6.49751e-07    44             1.21569e+02    
  (3, 1)      1.65842e-04    1.51087e-07    20             2.72393e+02    
  (4, 0)      1.28085e-04    1.58796e-08    11             1.69443e+02    
--------------------------------------------------------------------------------
Samples will be updated according to
-----------------------------
  index       N              
-----------------------------
  (0, 0)      16289
  (0, 1)      973
  (0, 2)      204
  (0, 3)      64
  (0, 4)      7
  (1, 0)      979
  (1, 1)      247
  (1, 2)      84
  (1, 3)      21
  (2, 0)      252
  (2, 1)      94
  (2, 2)      44
  (3, 0)      44
  (3, 1)      14
  (4, 0)      6
-----------------------------
Taking 9 additional samples at index (0, 3)...
Taking 32 additional samples at index (1, 2)...
Taking 1 additional sample at index (1, 3)...
Taking 1 additional sample at index (2, 1)...
Taking 29 additional samples at index (2, 2)...
Shape of the index set:
  ◼ 
  ◼ ◼ 
  ◼ ◼ ◼ 
  ◼ ◼ ◼ ◼ 
  ◼ ◼ ◼ ◼ ◼ 
--------------------------------------------------------------------------------
  index       E              V              N              W              
--------------------------------------------------------------------------------
  (0, 0)      7.71469e-02    6.01058e-03    25520          8.00000e+00    
  (0, 1)      5.45626e-03    1.24835e-04    1552           4.66274e+01    
  (0, 2)      1.59881e-03    9.36213e-06    294            8.01966e+01    
  (0, 3)      4.24411e-04    1.25601e-06    64             1.21569e+02    
  (0, 4)      9.18491e-05    1.89814e-08    16             1.69443e+02    
  (1, 0)      4.69202e-03    1.26412e-04    1786           4.66274e+01    
  (1, 1)      4.01388e-04    1.72171e-05    342            9.98823e+01    
  (1, 2)      1.72962e-04    3.20079e-06    84             1.89648e+02    
  (1, 3)      2.47017e-05    3.23237e-07    21             2.72393e+02    
  (2, 0)      1.01867e-03    1.43615e-05    407            8.01966e+01    
  (2, 1)      1.68391e-04    4.67795e-06    94             1.89648e+02    
  (2, 2)     -1.51018e-04    1.43917e-06    44             2.15217e+02    
  (3, 0)      4.38346e-04    6.49751e-07    44             1.21569e+02    
  (3, 1)      1.65842e-04    1.51087e-07    20             2.72393e+02    
  (4, 0)      1.28085e-04    1.58796e-08    11             1.69443e+02    
--------------------------------------------------------------------------------
Checking convergence...
  ==> Rates: α ≈ ( 1.680, 1.959), β ≈ ( 4.334, 4.095), γ ≈ ( 0.618, 0.618).
  ==> Variance of the estimator ≈ 6.84759e-07.
  ==> Bias of the estimator ≈ 5.59077e-05.
  ==> Non-trivial MSE splitting parameter ≈ 0.99.
Shape of the index set:
  ◼ 
  ◼ ◼ 
  ◼ ◼ ◼ 
  ◼ ◼ ◼ ◼ 
  ◼ ◼ ◼ ◼ ◼ 
--------------------------------------------------------------------------------
  index       E              V              N              W              
--------------------------------------------------------------------------------
  (0, 0)      7.71469e-02    6.01058e-03    25520          8.00000e+00    
  (0, 1)      5.45626e-03    1.24835e-04    1552           4.66274e+01    
  (0, 2)      1.59881e-03    9.36213e-06    294            8.01966e+01    
  (0, 3)      4.24411e-04    1.25601e-06    64             1.21569e+02    
  (0, 4)      9.18491e-05    1.89814e-08    16             1.69443e+02    
  (1, 0)      4.69202e-03    1.26412e-04    1786           4.66274e+01    
  (1, 1)      4.01388e-04    1.72171e-05    342            9.98823e+01    
  (1, 2)      1.72962e-04    3.20079e-06    84             1.89648e+02    
  (1, 3)      2.47017e-05    3.23237e-07    21             2.72393e+02    
  (2, 0)      1.01867e-03    1.43615e-05    407            8.01966e+01    
  (2, 1)      1.68391e-04    4.67795e-06    94             1.89648e+02    
  (2, 2)     -1.51018e-04    1.43917e-06    44             2.15217e+02    
  (3, 0)      4.38346e-04    6.49751e-07    44             1.21569e+02    
  (3, 1)      1.65842e-04    1.51087e-07    20             2.72393e+02    
  (4, 0)      1.28085e-04    1.58796e-08    11             1.69443e+02    
--------------------------------------------------------------------------------
Convergence reached. RMSE ≈ 8.29388e-04.
--------------------------------------------------------------------------------
*** MultilevelEstimators.jl @2018-07-09T15:31:40.049
*** Successfull termination
--------------------------------------------------------------------------------
```

### Multilevel Quasi-Monte Carlo Methods

When using plain Monte Carlo (MC) sampling, the optimal complexity rate of multilevel methods is <img src="https://latex.codecogs.com/svg.latex?\Large&space;1/\sqrt{N}"/>.  Using MIMC, we are able to achieve this optimal convergence rate for a large class of problems. However, even this optimal complexity rate is in practice much too expensive. Therefore, we couple our estimators with so-called Quasi-Monte Carlo (QMC) point set generators, that can achieve an integration error like <img src="https://latex.codecogs.com/svg.latex?\Large&space;1/N^\alpha"/>, with <img src="https://latex.codecogs.com/svg.latex?\Large&space;\alpha"/> > 0.5. In particular, we will focus on randomly shifted rank-1 lattice rules as implemented in the [QMC.jl](https://github.com/PieterjanRobbe/QMC.jl)-module.

To obtain an error bound on our estimator, we cannot simply take the variance over all samples, since the points are chosen deterministically. Instead, we must set up multiple copies of the same lattice rule, where each point set is shifted by a uniformly distributed random shift. The variance over these uncorrelated estimators can be used to construct a confidence interval the usual way.

Similarly, we can simulate the PDE with random coefficients from the previous section as

```julia
srand(2018) # for reproducability
estimator = init_lognormal_diffusion_mlqmc()
h = run(estimator,1e-3)
```

The module can also work with custom defined point generators, for example, the higher-order polynomial lattice rules from [QMC.jl](https://github.com/PieterjanRobbe/QMC.jl).

### Other Methods

It is also possible to combine MIMC with QMC, see [1]. The latest additons to this package are Adaptive Multi-Index Monte Carlo [2] and Multigrid Multilevel Monte Carlo [3].

### Parallel Sampling

Of course, as for any Monte Carlo method, all samples can be taken in parallel. The most convenient way is wrapping your code into a custom module, see, e.g., the [applications](https://github.com/PieterjanRobbe/MultilevelEstimators.jl/tree/master/applications/SPDE) folder.

Start julia in parallel or execute `addprocs(p)`, and include your module.

## Detailed description

A full list of parameters that can be specified is given below
<center>

| Parameter name         | Type                                      | Description                                    |
| ---------------------- | ----------------------------------------- | ---------------------------------------------- |
| `method` | <div>`SL`/`ML`</div><div>`TD`/`FT`/`HC`/`AD`</div> | Required. Index set to be used. Can be single level, multilevel, total degree, full tensor, hyperbolic cross, or adaptive. |
| `number_generator` |<div> `Uniform`/`Normal`+</div><div>`MC`/`QMC`+`Generator`</div> | Required. Point generator to be used. Can be Uniform (between `lb` and `ub`) or Normal, and Monte Carlo or Quasi-Monte Carlo|
| `sample_function`       | `Function`                                | Required. Function that specifies the quantity of interest or QOI. Inputs are an `Index` or `Level` to sample on, a `Vector` of random numbers, and an optional object with user data. This last input is usefull to pass objects (such as a `GaussianRandomField` or a Finite Element mesh) to the sample function. The mean of this QOI will be estimated. |
| `nb_of_warm_up_samples` | `Int64` | Number of warm-up samples. Default is `20` for Monte Carlo sampling, and `1` for Quasi-Monte Carlo. |
| `nb_of_qoi` | `Int64 ` | Number of QOI's that will be returned by the `sample_function`. Default is `1`. |
| `continuate` | `Bool` | Boolean. When equal to `true`, we run the simulation for a sequence of larger tolerances to obtain better estimates for the model parameters. Default is `false`.|
| `ntols` | `Int64` | Number of tolerances to be used in the continuation algorithm. Only used when `continuate` is `true`. Default is `10`. |
| `p0` | `Float64` | Parameter in the continuation algorithm. Determines the ratio between the different larger tolerances in the continuation. Default is `1.5`. |
| `name` | `String` | Problem name to be used when saving the simulation results. Default is `""`.|
| `folder` | `String` | Folder where the simulation results will be saved. Default is the current folder (returned by `cd()`)|
| `store_samples` | `Bool` | Option to store all samples taken. Usefull when plotting a probabiltiy distribution function of the QOI. Default is `false`.|
| `user_data` | `Any` | Option to add a user-defined object to the `Estimator`. Will be passed to the sample function. Default is `nothing`.|
| `verbose` | `Bool` | Log information to the terminal when `true`. Default is `false`.|
| `cost_model` | `Function` | A function that defines the cost to take a sample at a certain `level` or `index`. A predefined geometric cost model is provided as `(index) -> geometric_cost_model(index, M, c)`. When no cost model is provided, the actual run time is used. Beware of Julia's LLVM compiler though!|
| `conservative_bias_estimate` | `Bool` | When set to `false`, we use only the two finest levels to perform a regression and compute the bias. This is closer to the original method of Mike Giles. However, when the estimates of the expected value is poor, this might result in an inaccurate bias estimate and an early termination. |
| `max_level` | `Int64` | Maximum parameter `L` to be used in index set definition. The simulation will terminate whenever this value is reached, and the estimate so far will be returned. Default is `100`. |
| `max_search_space` | `IndexSet` | Shape of the search space in the adaptive MIMC algorithm, evaluated with `max_level`. Default is `TD`. |
| `do_regression` | `Bool` | When set to `true`, we compute regression on the number of samples on the finer levels instead of taking `nb_of_warm_up_samples` samples. Default is `true`. |
| `do_splitting` | `Bool` | When set to `true`, we compute a nontrivial splitting of the MSE. Default is `true`. |
| `parallel_sample_function` | `Function` | Usefull when defining your own parallel sample function. |
| `sample_multiplication_factor` | `Float64` | Determines the multiplication factor in the QMC *doubling* algorithm. Default is `2`. When `sample_multiplication_factor <=1 `, we compute one extra point at a time. |

</center>

[1] Robbe, P., Nuyens, D., Vandewalle, S. *A Multi-Index Quasi-Monte Carlo Algorithm for Lognormal Diffusion Problems.* SIAM Journal on Scientific Computing 39.5 (2017): S851-S872.

[2] Robbe, P., Nuyens, D., Vandewalle, S. *A Dimension-Adaptive Multi-Index Monte Carlo Method Applied to a Model of a Heat Exchanger.* Proceedings of the Twelfth International Conference on Monte Carlo and Quasi-Monte Carlo Methods in Scientific Computing (MCQMC2016).

[3] Robbe, P., Nuyens, D., Vandewalle, S. *Recycling Samples in the Multigrid Multilevel (Quasi-) Monte Carlo Method.* ArXiv preprint [arXiv:1806.05619](https://arxiv.org/pdf/1806.05619.pdf) (2018).
