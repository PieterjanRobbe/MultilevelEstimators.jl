# Multilevel Monte Carlo for Julia
[![Build Status](https://travis-ci.org/PieterjanRobbe/MultilevelEstimators.jl.png)](https://travis-ci.org/PieterjanRobbe/MultilevelEstimators.jl)

This module provides a versatile implementation of Multilevel Monte Carlo (MLMC) methods and Multi-Index Monte Carlo (MIMC) methods, and their Quasi-Monte Carlo (QMC) counterparts, Multilevel Quasi-Monte Carlo (MLQMC) and Multi-Index Quasi-Monte Carlo (MIQMC). The module is mainly aimed at solving partial differential equations (PDEs) with random coefficients, but can be used for any forward uncertainty quantification (UQ) problem by specifying the appropriate sample function.

The module has been tested for challenging 3D problems with up to 5000 uncertainties, and the parallel implementation was tested on a compute node with 24 cores.

Recent additions are Continuation Multilevel Monte Carlo (CMLMC) and Adaptive Multi-Index Monte Carlo (AMIMC).

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

In most mathematical models, parameters or coefficients are unknown or subject to uncertainty, particularly due to lack of data or measurements. Often, these problems involve the computation of a *quantity of interest* as the expected value over the uncertain input parameters. The classical sample-based approach then chooses `N` realisations of the uncertain parameters and approximates this expected value as a sample average. The Multilevel Monte Carlo (MLMC) method improves the error versus work complexity rate of the classical approach by using models  with different levels of accuracy. These models are called *levels*. The main idea is to write the approximation to G at the most accurate level `L` as a telescoping sum

<p align="center">
<img src="https://github.com/PieterjanRobbe/MultilevelEstimators.jl/blob/master/figures/eq_telsum.png" width="650" align="middle">
</p>

Note that `l=0` is the least accurate model. Hence, instead of approximating the expected value of the quantity of interest on the finest mesh, the MLMC method approximates differences &Delta;G at different levels `l`. If the variance of these differences goes sufficiently fast to zero as `l` increases, most samples are taken at models with low accuracy, hence low cost. Typically, only very few samples are needed at the finest mesh.


### A First Simple MLMC Example

PDEs with random coefficients typically describe cell movements, fluid pressures, temperatures or flow rates, all subject to uncertainty. For example, a model for fluid flow through porous media, commonly used in geophysics, is the elliptic PDE

<p align="center">
<img align="middle" src="https://github.com/PieterjanRobbe/MultilevelEstimators.jl/blob/master/figures/eq_PDE.png" width="250">
</p>

with random diffusion coefficient `k(x,ω)`. `k` represents the uncertain permeability of the porous material, and `p` is the unknown pressure head. The problem is physically meaningful if `k > 0`, hence, usually `k(x,ω) = exp(Z(x,ω))`, where `Z` is an underlying random field with known characteristics. For example, we may consider the Gaussian random field `Z(x,ω)` with Mat&eacute;rn kernel

<p align="center">
<img align="middle" src="https://github.com/PieterjanRobbe/MultilevelEstimators.jl/blob/master/figures/eq_matern.png" width="550">
</p>

Here, `p` denotes the usual p-norm. Let us choose

```julia
λ = 1. # correlation length
σ = 1. # variance
ν = 0.5 # smoothness
p = 1. # 1-norm
```
The Mat&eacute;rn kernel is predefined in the package, other kernels can be defined by specifying a covariance function.

```julia
myMaternKernel = MaternKernel(λ,σ,ν,p)
```
A typical sample of the lognormal random field is shown below:

<p align="center">
<img align="middle" src="https://github.com/PieterjanRobbe/MultilevelEstimators.jl/blob/master/figures/fig_matern.png" width="450">
</p>

When dealing with PDEs, the natural sequence of models with different accuracy, required in the Multilevel Monte Carlo method, can be obtained by varying the discretization mesh. For example, one can use a Finite Volume (FV) method and consider meshes with an increasing number of FV cells as the level `l` increases. We provide a simple implementation of the FV method using rectangular cells in

```julia
include("../2d/epde2.jl")
```

Let us consider the PDE in 2D on [0,1]&sup2;

```julia
d = 2 # physical dimension
```
and impose a geometrical grid hierarchy as follows:

<center>

| Level         | Number of FV cells           |
| ------------- |:-------------:|
| `0`           | `4`x`4` |
| `1`           | `8`x`8` |
| `2`           | `16`x`16` |
| `l`          | `(4*2^l)`x`(4*2^l)` |

</center>

Samples of a random field can be taken using the KL-expansion (see [1] for details). `MultilevelEstimators.jl` provides a KL expansion for random fields in arbitrary dimension with given covariance function `C`. We only need to specify the number of terms in the expansion as

```julia
s = 100 # 100 KL terms
```

Next, let us define a quantity of interest to compute. This function must take 3 input arguments: a vector of random variables used in the KL expansion, `xi`, an index (or level) `j` and a Settings-object that contains all simulation details specified above.

```julia
using Interpolations

function myQuantityOfInterest{T<:AbstractFloat,d,S<:Sampler}(xi::Vector{T},level::Index{d},sampler::S)
  kl = compose(sampler.gaussianFieldSampler,xi,level)::Matrix{T} # apply KL expansion
  m = 4.*2.^ level.indices # grid sizes
  mx, my = length(level.indices) < 2 ? tuple(repeat(m,inner=[2])...) : tuple(m...)
  k = exp(reshape(kl,(mx,my)))
  p = epde2(k) # solve the deterministic PDE (flow cell geometry)
  vx = 1/2/mx:1/mx:1-1/2/mx # FV computes solution in cell centers
  vy = 1/2/my:1/my:1-1/2/my
  myInterpolator = interpolate((vx,vy), p, Gridded(Linear()))
  return myInterpolator[0.5,0.5]
end
```
This code defines a quantity of interest `myQuantityOfInterest` that is a point evaluation of the pressure in the middle of the domain, `[0.5,0.5]`.

To set up a full simulation, we must provide the level hierarchy (or index set, see below), a number generator (MC or QMC, see below), a quantity of interest, here `myQuantityOfInterest`, and possibly a Gaussian field sampler (such as the KL expansion above).

```julia
myIndexSet = ML() # set-up for multilevel
myNumberGenerator = GaussianMCgenerator(s) # Monte Carlo sampler
myGaussianFieldSampler = KLExpansion(myMaternKernel,d,s,m0=4,maxL=6)
```

In the KL expansion, the random numbers are standard normal random variables. The function `KLExpansion` will precompute and store all `s` eigenfunctions on all `maxL+1` levels. We also provide a uniform random number generator, UniformMCgenerator, that generates points inside the hypercube `[0,1]^s`.

Next, create a `Dict` that holds all settings

```julia
myDict = Dict(
	"indexSet" => myIndexSet,
	"numberGenerator" => myNumberGenerator,
	"sampleFunction" => myQuantityOfInterest,
	"gaussianFieldSampler" => myGaussianFieldSampler
)
```
The `setup` function reads all settings from the `Dict` and returns an object of type `::Sampler`.

```julia
mySampler = setup(myDict)
```
Next, perform an MLMC simulation using the `simulate` function

```julia
srand(2016) # reset random number generator for reproducibility
TOL = 1e-3
t = simulate(mySampler,TOL)
```
The return value `t` is an array containing the elapsed time in seconds. The output of the simulation is

```
--------------------------------------------------------------------------------
*** MultilevelEstimators.jl @2016-11-08T12:00:00
*** Simulating myQuantityOfInterest
*** Using a Multilevel index set, not continuating 
*** absTOL = 1.000e-03 / relTOL = Inf (failure probability of 0.10)
--------------------------------------------------------------------------------
*** currently running at L = 0...
>> taking 32 samples at level 0...
   *** initial bias = 4.799021e-01 
   *** splitting = 5.000000e-01 
>> taking 180617 samples at level 0...
*** currently running at L = 1...
>> taking 32 samples at level 1...
   *** initial bias = 1.747457e-02 
   *** splitting = 5.000000e-01 
>> taking 53475 samples at level 1...
>> taking 77275 samples at level 0...
*** currently running at L = 2...
>> taking 32 samples at level 2...
   *** initial bias = 6.664800e-04 
   *** splitting = 5.000000e-01 
>> taking 258 samples at level 1...
>> taking 10147 samples at level 2...
>> taking 34275 samples at level 0...
*** error estimate is 5.045027e-04 + 2.270594e-04 = 7.315622e-04 
*** result is 4.995982e-01 ±(1.289122e-01, 2.578244e-01, 3.867367e-01) 
-------------------------------------------------------------------------------- 
 level       E              V               N               W                  
-------------------------------------------------------------------------------- 
  [0]         4.99769e-01    1.501859e-02    292199          1.000000e+00 
  [1]         5.65267e-05    1.437525e-03    53765           2.828427e+00 
  [2]        -2.27059e-04    1.622424e-04    10179           8.000000e+00
```

Hence, the expected value of the quantity of interest is 0.500 (up to `TOL=1e-3` accurate) and its 3σ-confidence interval is (0.113,0.886). Note that the exact expected value of this problem, `0.5`, is known exactly. The MLMC algorithm used three levels and took 292199 PDE solves at a grid with `16` cells, 53765 PDE solves at a grid with `64` cells, and 10179 PDE solves at a grid with `256` cells.

The current state of the Sampler is

```julia
mySampler
********************************************************************************
*                                SAMPLER STATUS                                *
********************************************************************************
indexSet        | Multilevel index set 
numberGenerator | 100-dimensional Gaussian Monte Carlo sampler with λ= 0.5 
sampleFunction  | myQuantityOfInterest 

...

------------------------------------------------------ 
  level           sample size           time/sample   
------------------------------------------------------
  [0]             (1,292199,1)            1.561920e-04
  [1]             (1,53765,1)             7.776886e-04
  [2]             (1,10179,1)             1.073885e-03
```

Note that there are different ways to call `simulate` on the `Sampler` type:

```julia
simulate(mySampler, absTOL)
simulate(mySampler, absTOL, failProb)
simulate(mySampler, absTOL=..., relTOL=..., failProb=...)
```

A full list of options that can be passed to the Sampler is given below. We also provide a `reset` function to clear the samples already taken:

```julia
reset(mySampler)
```

### Multi-Index Monte Carlo

Multi-Index Monte Carlo (MIMC) is a generalisation of MLMC where the levels are replaced by multi-indices. Whereas MLMC uses levels that refine in both `x`- and `y`-direction simultaneously, the MIMC method will allow for grids that refine only in `x` or only in `y`. For non-isotropic examples, where the MLMC method would perform badly, MIMC allows to again achieve the optimal convergence rate of error versus work. An illustration of these grids is shown below:

<p align="center">
<img align="middle" src="https://github.com/PieterjanRobbe/MultilevelEstimators.jl/blob/master/figures/fig_mimc.png" width="450">
</p>

The classical MLMC method would only consider meshes on the main diagonal, whereas MIMC considers a subset of all grids depicted. This subset is called the *index set*. As for the Smolyak construction in Sparse Grids, it can be shown that an index set consisting of the grids in the left upper triangle (the `TD` case) are, in some sense, optimal, see [2]. We provide 6 different index sets, `SL` (single level), `ML` (multilevel), `FT` (full tensor), `TD` (total degree), `HC` (hyperbolic cross) and `AD`. This last index set is an Adaptive Multi-Index Monte Carlo (AMIMC) method, see below.

Let us solve the model problem from the previous section using MIMC with TD index sets.

```julia
pd = 2 # physical dimension of the problem
d = 2 # dimension of the index set

indexset = TD(d) # setup for Multi-Index with TD index set

myDict = Dict(
    "indexSet" => myIndexSet,
    "numberGenerator" => myNumberGenerator,
    "sampleFunction" => myQuantityOfInterest,
    "gaussianFieldSampler" => myGaussianFieldSampler
)

mySampler = setup(myDict)
```
Note that the `d`-argument in the `TD` indexset is now equal to the physical dimension of the problem (`pd=2`), since we want our index set to refine in both `x`- and `y`-direction. When `d=1` (`ML()`), the standard MLMC method is used.

```julia
srand(2016) # reset random number generator for reproducibility
TOL = 5e-3
t = simulate(mySampler,TOL)
```

The result of the simulation is

```
--------------------------------------------------------------------------------
*** MultilevelEstimators.jl @2016-11-08T12:00:00
*** Simulating myQuantityOfInterest
*** Using a 2-dimensional index set of type TD, not continuating 
*** absTOL = 1.000e-03 / relTOL = Inf (failure probability of 0.10)
--------------------------------------------------------------------------------
*** currently running at L = 0...
>> taking 32 samples at index [0,0]...
   *** initial bias = 0.000000e+00 
   *** splitting = 1.000000e+00 
>> taking 45131 samples at index [0,0]...
*** currently running at L = 1...
>> taking 32 samples at index [0,1]...
>> taking 32 samples at index [1,0]...
   *** initial bias = 1.065292e-02 
   *** splitting = 5.000000e-01 
>> taking 17347 samples at index [0,1]...
>> taking 223687 samples at index [0,0]...
>> taking 45793 samples at index [1,0]...
*** currently running at L = 2...
>> taking 32 samples at index [2,0]...
>> taking 32 samples at index [0,2]...
>> taking 32 samples at index [1,1]...
   *** initial bias = 2.792704e-03 
   *** splitting = 5.000000e-01 
>> taking 10536 samples at index [2,0]...
>> taking 1493 samples at index [0,1]...
>> taking 58830 samples at index [0,0]...
>> taking 11116 samples at index [1,0]...
>> taking 2358 samples at index [0,2]...
>> taking 1691 samples at index [1,1]...
*** error estimate is 5.069647e-04 + 1.462006e-04 = 6.531653e-04 
*** result is 4.999281e-01 ±(1.289106e-01, 2.578212e-01, 3.867318e-01) 
-------------------------------------------------------------------------------- 
 index       E              V               N               W                  
--------------------------------------------------------------------------------
  [0,0]       5.00099e-01    1.503188e-02    327680          1.000000e+00 
  [0,1]      -8.29571e-05    1.414546e-04    18872           2.828427e+00 
  [1,0]      -4.80779e-05    1.283251e-03    56941           2.828427e+00 
  [0,2]      -1.51620e-05    5.751510e-06    2390            8.000000e+00 
  [1,1]      -7.78118e-05    4.032565e-06    1723            8.000000e+00 
  [2,0]       5.32269e-05    1.515711e-04    10568           8.000000e+00 
```

Hence, the expected value of the quantity of interest is again 0.500 (up to `TOL=1e-3` accurate) and its 3σ-confidence interval is (0.113,0.886). The MIMC algorithm uses six indices with 327680 PDE solves at a grid with `16` cell, 18872 + 56941 PDE solves at a grid with `32` cells 2390 + 1723 + 10568 PDE solves at a grid with `64` cells.

### Multilevel Quasi-Monte Carlo Methods

When using plain Monte Carlo (MC) sampling, the optimal complexity rate of multilevel methods is *error &sim; work*^(-&frac12;). Using MIMC, we are able to achieve this optimal convergence rate for a large class of problems. However, even this optimal complexity rate is in practice much too expensive. Therefore, we couple our estimators with so-called Quasi-Monte Carlo (QMC) point set generators, that can achieve an integration error like *error &sim; N^(-&alpha;)*, with α > 0.5. In particular, we will focus on randomly shifted rank-1 lattice rules as implemented in the [`QMC.jl`](https://github.com/PieterjanRobbe/QMC.jl)-module.

To obtain an error bound on our estimator, we cannot simply take the variance over all samples, since the points are chosen deterministically. Instead, we must set up multiple copies of the same lattice rule, where each point set is shifted by a uniformly distributed random shift. The variance over these `q` estimators can be used to construct a confidence interval the usual way. For example, let us use a lattice rule with `16` shifts

```julia
q = 16 # number of shifts
myNumberGenerator = GaussianQMCgenerator(s,q) # Quasi-Monte Carlo sampler
```
Note that also a uniform sampler, `UniformQMCgenerator`, is available.

To complete the simulation, execute

```julia
myDict = Dict(
    "indexSet" => myIndexSet,
    "numberGenerator" => myNumberGenerator,
    "sampleFunction" => myQuantityOfInterest,
    "gaussianFieldSampler" => myGaussianFieldSampler
)

mySampler = setup(myDict)

srand(2016) # reset random number generator for reproducability
TOL = 1e-3
t = simulate(mySampler,TOL)
```

This results in

```
--------------------------------------------------------------------------------
*** MultilevelEstimators.jl @2016-11-08T19:32:12
*** Simulating myQuantityOfInterest
*** Using a Multilevel index set, not continuating 
*** absTOL = 1.000e-03 / relTOL = Inf (failure probability of 0.10)
--------------------------------------------------------------------------------
*** currently running at L = 0...
>> taking 16x2 samples at level 0...
   *** initial bias = 4.899834e-01 
   *** splitting = 5.000000e-01 
>> taking 16x113 samples at level 0...
*** currently running at L = 1...
>> taking 16x2 samples at level 1...
   *** initial bias = 1.974494e-03 
   *** splitting = 5.000000e-01 
>> taking 16x37 samples at level 1...
>> taking 16x19 samples at level 0...
>> taking 16x25 samples at level 1...
*** currently running at L = 2...
>> taking 16x2 samples at level 2...
   *** initial bias = 2.479182e-03 
   *** splitting = 5.000000e-01 
>> taking 16x18 samples at level 2...
>> taking 16x33 samples at level 0...
>> taking 16x64 samples at level 1...
>> taking 16x12 samples at level 2...
*** error estimate is 5.861228e-04 + 1.626384e-04 = 7.487612e-04 
*** result is 5.009067e-01 ±(1.259022e-01, 2.518044e-01, 3.777066e-01) 
--------------------------------------------------------------------------------
  level       E              V               N               W                 
-------------------------------------------------------------------------------- 
  [0]         5.00254e-01    1.424236e-02    2672            1.000000e+00 
  [1]         4.90160e-04    1.447303e-03    2048            2.828427e+00 
  [2]         1.62638e-04    1.617048e-04    512             8.000000e+00 
```

The MLQMC method only needs 2672 samples of the PDE on a `16`-point mesh, 2048 samples on a `64`-point mesh, and 512 on the finest mesh. Compare this to the 292199  samples at the coarsest mesh for MLMC. We can check that indeed the sum of the statistical error and bias is smaller than `TOL = 1e-3`.

The module can also work with custom defined point generators, for example, the higher-order polynomial lattice rules from [QMC.jl](https://github.com/PieterjanRobbe/QMC.jl).

### Other Methods

The latest additons to the library are Adaptive Multi-Index Monte Carlo and Continuation Multilevel Monte Carlo. In Adaptive MIMC, the index set is built adaptively. This has the adavantage that one does not need insight into the structure of the problem across all indices. The index sets from AMIMC can be very different from the known FT, TD or HC type. 

In Continuation MLMC simulation, the problem parameters are estimated on the fly using a fit through the available statistics from the ensemble. This way, one avoids taking a lot of warm-up samples on the fine grids, where samples are expensive. Continuation also works with general index sets, although theory is lacking.

### Multiple Quantities of Interest

It is also possible to compute multiple quantities of interest in the same simulation. For example, if instead of a single point evaluation of the pressure, we would like to find the statistics of the pressure in `121` points `(0:0.1:1,0:0.1:1)`. Therefore, we must define the `Dict` as

```julia
myDict = Dict(
	"indexSet" => myIndexSet,
	"numberGenerator" => myNumberGenerator,
	"sampleFunction" => myQuantityOfInterest,
	"gaussianFieldSampler" => myGaussianFieldSampler,
	"Z"=>121
)
```
and adapt the quantity of interest to interpolate in multiple points:

```julia
function myQuantityOfInterest{T<:AbstractFloat,d,S<:Sampler}(xi::Vector{T},level::Index{d},sampler::S)
  kl = compose(sampler.gaussianFieldSampler,xi,index)::Matrix{T} # apply KL expansion
  m = [4,4].*2.^ level.indices # grid sizes
  mx, my = length(level.indices) < 2 ? tuple(repeat(m,inner=[2])...) : tuple(m...)
  k = exp(reshape(kl,(mx,my)))
  p = epde2d(k) # solve the deterministic PDE (flow cell geometry)
  vx = 1/2/mx:1/mx:1-1/2/mx # FV computes solution in cell centers
  vy = 1/2/my:1/my:1-1/2/my
  myInterpolator = interpolate((vx,vy), p, Gridded(Linear()))
  z = collect(0:0.1:1)
  return myInterpolator[z,z]
end
```

The result of the simulation will contain the expected value, variance and error estimate in each interpolation point.

### Parallel Sampling

Of course, as for any Monte Carlo method, all samples can be taken in parallel. The most convenient way is wrapping your code into a custom module, see the [test](https://github.com/PieterjanRobbe/MultilevelEstimators.jl/tree/master/test) folder.

```julia
module TestModule

	using Interpolations
	
	using Reexport
    @reexport using MultilevelEstimators
    
	export myQuantityOfInterest

	include("path/to/solvers/epde2.jl")
	
end
```

Start julia in parallel or execute `addprocs(p)`, and include your module:

```julia
push!(LOAD_PATH,".") # make sure the test module can be found

using TestModule # include the test module
```

Samples will be taken in parallel now.

## Detailed description

A full list of parameters that can be specified is given below
<center>

| Parameter name         | Type                                      | Description                                    |
| ---------------------- | ----------------------------------------- | ---------------------------------------------- |
| `indexSet`             | <div>`SL`/`ML`</div><div>`TD`/`FT`/`HC`/`AD`</div> | Index set to be used. Can be single level, multilevel, total degree, full tensor, hyperbolic cross, or adaptive. |
| `numberGenerator`      |<div> `Uniform`/`Gaussian`+</div><div>`MC`/`QMC`+`Generator`</div> | Point generator to be used. Can be Uniform (between `lb` and `ub`) or Gaussian, Monte Carlo or QuasiMonte Carlo|
| `sampleFunction`       | `Function`                                | Function that specifies the quantity of interest or QOI. Inputs are a `Vector` of random numbers, an `Index` and a `Sampler`. The mean of this QOI will be estimated. |
| `maxL`                 | `Int64`                                   | Maximum parameter L to be used in index set definition. The simulation will terminate whenever this value is reached, and return the result achieved so far. Default is `6`.|
| `γ`                    | `Vector{Float64}`                         | Parameters in the cost model `cost(index)~prod(2^(γ*index))`. Default is `1.5*ones(Float64,d)`. |
| `Z`                    | `Int64 `                                  | Number of QOIs that will be returned by the `sampleFunction`. Default is `1`. |
| `Nstar`                | `Int64 `                                  | Number of warm-up samples. Default is `32/nshifts(numberGenerator)`, so, `32` for plain Monte Carlo sampling. |
| `gaussianFieldSampler` | `KLExpansion`                             | Usefull when using random fields in QOI definition. Can be something of type `GaussianFieldSampler`, or `Vector{GaussianFieldSampler}`. Default is `EmptySampler()`. |
| `useTime`              | `Bool`                                    | Boolean. Is equal to `true` when the estimated actual simulation time is used in the computation of the number of samples. Use this if no cost model is available. Default is `false`.|
| `safety`               | `Bool`                                    | Boolean. Is equal to `true` when the variance of the estimator is guaranteed to be smaller then `TOL^2/2`. `safety` must always be `true` when using QMC point generators. Default is `true`. |
| `continuate`           | `Bool`                                    | Boolean. Is equal to `true` when doing continuation. Default is `false`.|
| `nTOL`                 | `Int64`                                   | Number of tolerances to be used in the continuation algorithm. Only specify `nTOL` when `continuate` is `true`. Default is `10`. |
| `k`                    | `Tuple{Float64,Float64}`                  | Trust parameters for variance estimation in the CMLMC algorithm. Default is `(0.1,0.1)`. |
| `showInfo`             | `Bool`                                    | Boolean. Is equal to `true` when info must be printed to `ioStream`. Default is `true`. |
| `ioStream`             | `IO`                                      | Where to write output of the simualtion to. Default is `STDOUT`. |
| `storeSamples`         | `Bool`                                    | Boolean. Is equal to `true` when not only differences, but also the values at each index must be stored. Usefull for diagnostics. Default is `false`. |
| `procMap`              | `Dict{Index,Int64}`                       | Allows to specify the distribution of processors across all levels. Each entry in the Dict is an `Index` that specifies the number of processors to be used when sampling at that `Index`. Default is a `Dict` with all entries equal to `nprocs()` - the maximum number of processors available. | 
| `userType`             | `T`                                       | Allows to append a user-defined type to the `Sampler`. Usefull when for example, meshes of a FE approximation must be loaded. Default is `Void`. |

</center>

[1] Le Maître, O., Knio., O. M. *Spectral methods for uncertainty quantification*. Springer Science & Business Media, 2010.

[2] Haji-Ali, A.-L., Nobile, F., Tempone, R. *Multi-index Monte Carlo: when sparsity meets sampling*. Numerische Mathematik 132.4 (2016): 767-806.
