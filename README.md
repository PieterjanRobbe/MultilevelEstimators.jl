# Multilevel Monte Carlo for Julia

This module provides a versatile implementation of Monte Carlo (MC) methods, Multilevel Monte Carlo (MLMC) methods and Multi-Index Monte Carlo (MIMC) methods, and their Quasi-Monte Carlo (QMC) counterparts, Multilevel Quasi-Monte Carlo (MLQMC) and Multi-Index Quasi-Monte Carlo (MIQMC). The module is mainly aimed at solving partial differential equations (PDEs) with random coefficients, but can be used for any forward uncertainty quantification (UQ) problem by specifying the appropriate sample function.

The module has been tested for challenging 3D problems with up to 5000 uncertainties, and the parallel implementation was tested on a compute node with 24 cores.

## Installation

```julia
Pkg.clone("https://github.com/PieterjanRobbe/MultilevelEstimators.jl")
```

## Usage

### Multilevel Monte Carlo

In most mathematical models, parameters or coefficients are unknown or subject to uncertainty, particularly due to lack of data or measurements. Often, these problems involve the computation of a *quantity of interest* as the expected value over the uncertain input parameters. The classical sample-based approach then chooses `N` realisations of the uncertain parameters and approximates this expected value as a sample average. The Multilevel Monte Carlo (MLMC) method improves the error versus work complexity rate of the classical approach by using models  with different levels of accuracy. These models are called *levels*. Note that `l=0` is the least accurate model. The main idea is to write the approximation to G at the most accurate level `L` as a telescoping sum

<img src="https://github.com/PieterjanRobbe/MultilevelEstimators.jl/blob/master/figures/eq_telsum.png" width="650" align="middle">

.. image:: https://github.com/PieterjanRobbe/MultilevelEstimators.jl/blob/master/figures/eq_telsum.png
    :alt: telsum
    :width: 650
    :align: center

Hence, instead of approximating the expected value of the quantity of interest on the finest mesh, the MLMC method approximates differences &Delta;G at different levels `l`. If the variance of these differences goes sufficiently fast to zero as `l` increases, most samples are taken at models with low accuracy, hence low cost. Typically, only very few samples are needed at the finest mesh.


### A First Simple MLMC Example

PDEs with random coefficients typically describe cell movements, fluid pressures, temperatures or flow rates, all subject to uncertainty. For example, a model for fluid flow through porous media, commonly used in geophysics, is the elliptic PDE
<img align="middle" src="https://github.com/PieterjanRobbe/MultilevelEstimators.jl/blob/master/figures/eq_PDE.png" width="250">

with random diffusion coefficient `k(x,ω)`. `k` represents the uncertain permeability of the porous material, and `p` is the unknown pressure head. The problem is physically meaningful if `k > 0`, hence, usually `k(x,ω) = exp(Z(x,ω))`, where `Z` is an underlying random field with known characteristics. For example, we may consider the Gaussian random field `Z(x,ω)` with Mat&eacute;rn kernel

<img align="middle" src="https://github.com/PieterjanRobbe/MultilevelEstimators.jl/blob/master/figures/eq_matern.png" width="550">

Here, `p` denotes the usual p-norm. Let us choose

```julia
λ = 1. # correlation length
σ = 1. # variance
ν = 0.5 # smoothness
p = 1 # 1-norm
```
We can define the covariance function `C` using the `matern`-function:

```julia
C(x,y) = matern(λ,σ,ν,x,y)
```
A typical sample of the lognormal random field is shown below:
<img align="middle" src="https://github.com/PieterjanRobbe/MultilevelEstimators.jl/blob/master/figures/fig_matern.png" width="450">

When dealing with PDEs, the natural sequence of models with different accuracy, required in the Multilevel Monte Carlo method, can be obtained by varying the discretization mesh. For example, one can use a Finite Volume (FV) method and consider meshes with an increasing number of FV cells as the level `l` increases. We provide a simple implementation of the FV method using rectangular cells in

```julia
include("../2d/epde2.jl")
```

Let us consider the PDE in 2D on [0,1]&sup2;

```julia
d = 2 # physical dimension
```
and impose a geometrical grid hierarchy as follows:

| Level         | Number of FV cells           |
| ------------- |:-------------:|
| `0`           | `4`x`4` |
| `1`           | `8`x`8` |
| `2`           | `16`x`16` |
| `l`          | `(4*2^l)`x`(4*2^l)` |


Samples of a random field can be taken using the KL-expansion (see [1] for details). `MultilevelEstimators.jl` provides a KL expansion for random fields in arbitrary dimension with given covariance function `C`. We only need to specify the number of terms in the expansion as

```julia
s = 100 # 100 KL terms
```

Next, let us define a quantity of interest to compute. This function must take 3 input arguments: a vector of random variables used in the KL expansion, `xi`, an index (or level) `j` and a Settings-object that contains all simulation details specified above.

```julia
using Interpolations

function myQuantityOfInterest{d,T}(xi::Vector{T}, j::Index{d}, set::Settings)
	k = compose(set.gaussianFieldSampler,xi,j) # apply KL expansion
	m = set.m0.*2.^j.indices # grid sizes
	mx, my = length(j.indices) < 2 ? tuple(repeat(m,inner=[2])...) : tuple(m...)
	k = exp(reshape(k,(mx,my)))
	p = epde2(k)
	vx = 1/2/mx:1/mx:1-1/2/mx # FV computes solution in cell centers
	vy = 1/2/my:1/my:1-1/2/my
	myInterpolator = interpolate((vx,vy), p, Gridded(Linear()))
	return myInterpolator[0.5,0.5]
end
```
This code defines a quantity of interest `myQuantityOfInterest` that is a point evaluation of the pressure in the middle of the domain, `[0.5,0.5]`.

To set up a full simulation, we must provide the level hierarchy (or index set, see below), a number generator (MC or QMC, see below), a Gaussian field sampler (such as the KL expansion above) and the quantity of interest `myQuantityOfInterest`.

```julia
indexset = createIndexSet(ML,1) # setup for Multilevel
numberGenerator = GaussianMCgenerator(1,s) # Monte Carlo sampler
gaussianFieldSampler = createKLexpansion(d,λ,σ,ν,s,cov=C) # KLE
sampleFunction = myQuantityOfInterest 
```

In the KL expansion, the random numbers are standard normal random variables. We also provide a uniform random number generator, UniformMCgenerator, that generates points inside the hypercube `[0,1]^s`.

Next, create a `Settings`-object and a `Sampler`:

```julia
mySettings = Settings(indexset, numberGenerator, gaussianFieldSampler, sampleFunction)
mySampler = createSampler(1,mySettings)
```

and perform a simulation as

```julia
srand(2016) # reset random number generator for reproducibility
TOL = 1e-3
(E, V, A, B, splitting, S) = simulate(mySampler,TOL)
```
Here, `TOL` is the desired tolerance on the *root mean square error* (RMSE). The `simulate` function has 6 different outputs:

* `E`: the expected value of the quantity of interest, computed up to `TOL`
* `V`: the variance of the quantity of interest
* `A`: the estimated statistical error
* `B`: the estimated bias (the total error is then `A+B < TOL`)
* `splitting`: the error splitting between statistical error and bias, i.e., `A = splitting*TOL` and `B=(1-splitting)*TOL`
* `S`: the optimal number of samples to take at each level (a `Dict()`-object)

The result of the simulation is

```
*** currently running at K = 0...
>> Taking 16 samples at 1-dimensional index with values [0]...
   >>> B = 0.5097936451161305
   >>> splitting = 0.5
>> Taking 137560 samples at 1-dimensional index with values [0]...
>> Taking 124568 samples at 1-dimensional index with values [0]...
*** currently running at K = 1...
>> Taking 16 samples at 1-dimensional index with values [1]...
   >>> B = 0.0009228244808867983
   >>> splitting = 0.5
>> Taking 63049 samples at 1-dimensional index with values [1]...
>> Taking 5811 samples at 1-dimensional index with values [0]...
*** error estimate is 0.0004920494062513982 + 2.4416544003517002e-5 = 0.0005164659502549153
*** result is 0.5000289938207338 with a variance of 0.017320835721944317
```

Hence, the expected value of the quantity of interest is 0.500 (up to `TOL=1e-3` accurate) and its variance is 0.017. Note that the exact expected value of this problem is known analytically:`E=0.5`. The MLMC algorithm used only two levels and took 267955 PDE solves at a grid with `16` cells and 63065 PDE solves at a grid with `64` cells.

We also provide a `reset!` function to clear the samples already taken:

```julia
reset!(mySampler)
```

### Multi-Index Monte Carlo

Multi-Index Monte Carlo (MIMC) is a generalisation of MLMC where the levels are replaced by multi-indices. Whereas MLMC uses levels that refine in both `x`- and `y`-direction simultaneously, the MIMC method will allow for grids that refine only in `x` or only in `y`. For non-isotropic examples, where the MLMC method would perform badly, MIMC allows to again achieve the optimal convergence rate of error versus work. An illustration of these grids is shown below:

<img align="middle" src="https://github.com/PieterjanRobbe/MultilevelEstimators.jl/blob/master/figures/fig_mimc.png" width="450">

The classical MLMC method would only consider meshes on the main diagonal, whereas MIMC considers a subset of all grids depicted. This subset is called the *index set*. As for the Smolyak construction in Sparse Grids, it can be shown that an index set consisting of the grids in the left upper triangle (the `TD` case) are, in some sense, optimal, see [2]. We provide 5 different index sets, `ML` (multilevel), `FT` (full tensor), `TD` (total degree), `HC` (hyperbolic cross) and `AD`. This last index set is an Adaptive Multi-Index Monte Carlo (AMIMC) method.

Let us solve the model problem from the previous section using MIMC with TD index sets.

```julia
pd = 2 # physical dimension of the problem
d = 2 # dimension of the index set

indexset = createIndexSet(TD,d) # setup for Multilevel
numberGenerator = GaussianMCgenerator(d,s) # Monte Carlo sampler
gaussianFieldSampler = createKLexpansion(pd,λ,σ,ν,s,cov=C) # KLE

mySettings = Settings(indexset, numberGenerator, gaussianFieldSampler, sampleFunction)
mySampler = createSampler(d,mySettings)
```
Note that the `d`-argument in the setup functions is now equal to the physical dimension of the problem (`pd=2`), since we want our index set to refine in both `x`- and `y`-direction. When `d=1`, the standard MLMC method is used.

```julia
srand(2016) # reset random number generator
(E, V, A, B, splitting, S) = simulate(mySampler,TOL)
```

The result of the simulation is

```
*** currently running at K = 0...
>> Taking 16 samples at 2-dimensional index with values [0,0]...
   >>> B = 0.0
   >>> splitting = 1.0
>> Taking 34378 samples at 2-dimensional index with values [0,0]...
>> Taking 31142 samples at 2-dimensional index with values [0,0]...
*** currently running at K = 1...
>> Taking 16 samples at 2-dimensional index with values [0,1]...
>> Taking 16 samples at 2-dimensional index with values [1,0]...
   >>> B = 0.005631158516032551
   >>> splitting = 0.5
>> Taking 16748 samples at 2-dimensional index with values [0,1]...
>> Taking 223827 samples at 2-dimensional index with values [0,0]...
>> Taking 66146 samples at 2-dimensional index with values [1,0]...
*** error estimate is 0.0004963130465833095 + 0.0002511187661258237 = 0.0007474318127091333
*** result is 0.500104729201501 with a variance of 0.0172919184738401
```

Hence, the expected value of the quantity of interest is again 0.500 (up to `TOL=1e-3` accurate) and its variance is 0.017. The MIMC algorithm uses three indices with 289363 PDE solves at a grid with `16` cells and 16764 + 66162 at a grid with `32` cells. 

### Multilevel Quasi-Monte Carlo Methods

When using plain Monte Carlo (MC) sampling, the optimal complexity rate of multilevel methods is *error &sim; work*^(-&frac12;). Using MIMC, we are able to achieve this optimal convergence rate for a large class of problems. However, even this optimal complexity rate is in practice much too expensive. Therefore, we couple our estimators with so-called Quasi-Monte Carlo (QMC) point set generators, that can achieve an integration error like *error &sim; N^(-&alpha;)*, with α > 0.5. In particular, we will focus on randomly shifted rank-1 lattice rules as implemented in the [`QMC.jl`](https://github.com/PieterjanRobbe/QMC.jl)-module.

To obtain an error bound on our estimator, we cannot simply take the variance over all samples, since the points are chosen deterministically. Instead, we must set up multiple copies of the same lattice rule, where each point set is shifted by a uniformly distributed random shift. The variance over these `q` estimators can be used to construct a confidence interval the usual way. For example, let us use a lattice rule with `16` shifts

```julia
q = 16 # number of shifts
numberGenerator = GaussianQMCgenerator(1,s,q) # Quasi-Monte Carlo sampler
```
The `1` here indicates that we will use Multilevel Quasi-Monte Carlo (MLQMC). Note that also a uniform sampler, `UniformQMCgenerator`, is available.

To complete the simulation, execute

```julia
indexset = createIndexSet(ML,1) # setup for Multilevel
gaussianFieldSampler = createKLexpansion(d,λ,σ,ν,s,cov=C) # KLE
sampleFunction = myQuantityOfInterest

mySettings = Settings(indexset, numberGenerator, gaussianFieldSampler, sampleFunction)
mySampler = createSampler(1,mySettings)

srand(2016) # reset random number generator for reproducability
(E, V, A, B, splitting, S) = simulate(mySampler,TOL)
```

This results in

```
*** currently running at K = 0...
>> Taking 16x16 samples at 1-dimensional index with values [0]...
   >>> B = 0.4983587937565641
   >>> splitting = 0.5
>> Taking 16x85 samples at 1-dimensional index with values [0]...
*** currently running at K = 1...
>> Taking 16x16 samples at 1-dimensional index with values [1]...
   >>> B = 0.0006964829951229075
   >>> splitting = 0.5
>> Taking 16x36 samples at 1-dimensional index with values [1]...
>> Taking 16x39 samples at 1-dimensional index with values [0]...
*** error estimate is 0.0001937976494168214 + 0.00046909667415580606 = 0.0006628943235726275
*** result is 0.49971563624899434 with a variance of 0.01685275419659606
```

The MLQMC method only needs `q`x 140 = 2240 samples of the PDE on a `16`-point mesh, and `q`x 52 = 832 samples on a `64`-point mesh. We can check that indeed the sum of the statistical error and bias `A + B ` is smaller than `TOL = 1e-3`.

The extension of the module to higher-order polynomial lattice rules is ongoing work.

### Multiple Quantities of Interest

It is also possible to compute multiple quantities of interest in the same simulation. For example, if instead of a single point evaluation of the pressure, we would like to find the statistics of the pressure in `121` points `(0:0.1:1,0:0.1:1)`. Therefore, we must define the `Settings` object as

```julia
mySettings = Settings(indexset, numberGenerator, gaussianFieldSampler, sampleFunction, Z=121)
```
and adapt the quantity of interest to interpolate in multiple points:

```julia
function myQuantityOfInterest{d,T}(xi::Vector{T},j::Index{d},set::Settings)
  k = compose(set.gaussianFieldSampler,xi,j) # KL expansion
  m = set.m0.*2.^j.indices # grid sizes
  mx, my = length(j.indices) < 2 ? tuple(repeat(m,inner=[2])...) : tuple(m...)
  k = exp(reshape(k,(mx,my))) # take exponential
  p = epde2(k) # solve the deterministic PDE (flow cell geometry)
  vx = 1/2/mx:1/mx:1-1/2/mx
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

[1] Le Maître, O., Knio., O. M. *Spectral methods for uncertainty quantification*. Springer Science & Business Media, 2010.

[2] Haji-Ali, A.-L., Nobile, F., Tempone, R. *Multi-index Monte Carlo: when sparsity meets sampling*. Numerische Mathematik 132.4 (2016): 767-806.
