## estimator.jl : main Estimator type
#
# Representation of the main Estimator type.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for
# Multilevel Monte Carlo Methods

## Estimator ##
"""
    Estimator(index_set::AbstractIndexSet,
              sample_method::AbstractSampleMethod,
              sample_function::Function,
              distributions::Vector{<:Distribution}; kwargs...)

Create an `Estimator` with index set `index_set` and sample method `sample_method` to compute the expected value of the quantity of interest returned by `sample_function`, and where `distributions` is the uncertainty on the input parameters.

# Examples

See the [example section in the online documentation](@ref Example) for examples on how to use this function.

# General keywords

Different algorithmic options can be set by providing appropriate keyword arguments. A list of these keyword arguments and their description is provided below.

!!! note

    Not all keyword arguments can be used in combination with all `Estimator` types.

**`name`** -- a human readable name for the problem that must be solved. This name will be used when saving a [`History`](@ref) file with diagnostic information from each [`run`](@ref), and also when generating reports with the [`Reporter`](https://github.com/PieterjanRobbe/Reporter.jl)-package. Default is `UntitledEstimator`.

**`folder`** is where the [`History`](@ref) files with diagnostic information from every [`run`](@ref) will be saved. Default is the current directory.

**`save`** -- if set to `true`, MultilevelEstimators saves diagnostic information about the simulation in a separate file in the `.jld2` format. See [`JLD2`](https://github.com/JuliaIO/JLD2.jl) for more information about this file format. By default, this keyword is set to `true`.

**`verbose`** -- if set to `true`, MultilevelEstimators will print useful information about the [`run`](@ref) to `STDOUT`. Default is `true`.

**`nb_of_warm_up_samples`** is the number of warm-up samples to take at each level (or index) to obtain initial variance estimates. In combination with `do_regression = true` (default), this is only true for all levels < 2 (or indices with all coordinate directions < 2). On all other levels (or indices), this is the initial maximum number of regressed samples. When using [`MC`](@ref)-sampling, this value must be larger than or equal to 2, in order to be able to compute the variance. The default value is 20. When using [`QMC`](@ref)-sampling, this is the number of samples per random shift, and since we only consider the variance with repect to all random shifts, this value can be as small as 1 (default). Increase this number for more accurate initial variance estimation.   

**`nb_of_qoi`** is the number of quantities of interest returned by the `sample_function`. In this case, the simulation will be performed by controlling the root mean square error of the quantity of interest with the largest variance. The expected value of all other quantities of interest is guaranteed to be estimated more accurate or as accurate as the one with the largest variance. The default value is 1. The quantitiy of interest can optionally be provided using the key `qoi_with_max_var`, te ensure that the same quantity is used in every iteration. 

**`continuate`** -- if set to `true`, each call to [`run`](@ref) with a specified tolerance `tol` results in a simulation for a sequence of larger tolerances ``\\text{tol}_i = p^{(n-1-i)} \\text{tol}`` to get improved estimates for the variance and the bias. The default value is `true`.

**`continuation_mul_factor`** is the value ``p`` in the above formula. The default value is 1.2.

**`nb_of_tols`** is the value ``n`` in the above formula. The default value is 10.

**`save_samples`** -- if set to `true`, all samples of a [`run`](@ref) will be saved together with the diagnostic information. This is usefull for plotting histograms etc.

**`cost_model`** allows the user to provide a function that returns the cost of a sample on each index. For example, a geometric cost model could be provided with

```julia
γ = 1.5
Estimator(...; cost_model = level -> 2^(γ * level[1]))
```

When no cost model is provided (default), the actual run time on each level (or index) is used as a proxy for the cost model.

!!! tip

    Beware of precompilation! Using actual run times for the cost model might be very unreliable, especially on the coarsest level (or index), when the [`Estimator`](@ref) is [`run`](@ref) for the first time. In real production code, it is advised that the `run`-function is called once beforehand with a larger tolerance than the target accuracuy, to eliminate precompilation time.
    ```julia
    run(estimator, 10 * tol) # precompilation run
    run(estimator, tol)
    ```
 
**`nb_of_workers`** allows the user to specify the number of workers (processors) to be used for parallel sampling in a multicore environment. By default, this is equal to the value returned by `nworkers()`. `nb_of_workers` can also be given as function, specifying the number of workers on each level (or index). This is particularly usefull for expensive sample functions that require non-trivial load balancing.

**`nb_of_uncertainties`** is a function indicating the number of uncertainties to be used on each level (or each index). The default value is equal to the length of `distributions`, the vector of [`Distribution`](@ref Distribution)s specifying the input uncertainty of the parameters. This is useful for so-called level-dependent estimators, where a different number of parameters is activated on each level (or index). The length of `distributions` must be larger than or equal to the maximum value returned by this function.

**`max_index_set_param`** is the maximum index set size parameter that can be used in a simulation. Typically, a simulation will only be possible up to a finite index set size parameter, due to, for instance, memory constraints or prohibitively lare computational cost. For a multilevel simulation with [`ML`](@ref) this is the maximum number of levels. For a multi-index simulation with [`FT`](@ref), [`TD`](@ref), [`HC`](@ref) or [`ZC`](@ref), this is the maximum index set size parameter. For simulations with [`AD`](@ref) or [`U`](@ref), this is the maximum index set size to be used in the `max_search_space` (see below). Default is 10 times the number of dimensions of the index set. When this value is exceeded, the simulation will stop with a bias (and, consequently, a root mean square error) larger than the requested accuracy.

**`min_index_set_param`** is the minimum index set size parameter that should be used in a simulation. Useful for providing a minimum number of levels that should be used. This parameter should be at least 2 to allow for regression and estimation of the rates `alpha` and `beta. Default value is 2.

# Keywords valid for index set types [`SL`](@ref), [`ML`](@ref), [`FT`](@ref), [`TD`](@ref), [`HC`](@ref) and [`ZC`](@ref)

**`do_regression`** -- if set to `true`, the variance on the finer levels (or indices) will be guessed from the aready available variances on the coarser levels (or indices). These finer levels (or indices) must have a value larger than or equal to 2 in at least one of the coordinate directions. This is very powerfull when dealing with computationally demanding sample functions. Default: `true`.

**`do_mse_splitting`** -- if set to `true`, MultilevelEstimators will use a non-trivial splitting of the mean square error of the estimator: ``\\text{MSE} = \\theta \\text{var} + (1 - \\theta) \\text{bias}^2``. The value of the splitting parameter is very important to obtain an efficient estimator. Default is `true`.  

**`min_splitting`** is the minimum value for the mean square error splitting parameter. Default: 0.5.

**`max_splitting`** is the maximum value for the mean square error splitting parameter. Default: 0.99.

# Keywords valid for index set types [`AD`](@ref) and [`U`](@ref)

**`max_search_space`** determines, together with `max_index_set_param`, the collection of all indices that may be considered for further refinement. This is necessary for `AD`- and `U`-type index sets, that do not define such a set, by defintion. Default: TD.  

# Keywords valid for [`AD`](@ref) index sets

**`penalization`** is a factor added to the computation of the profit indicators in the adaptive algorithm. We observed that in some cases, the algorithm stretches its search space too much around the coordinate axes, especially when actual computation times are used (see `cost_model`). The penalization parameter 0<`p`<1 penalizes these search directions, in favor of more classically-shaped index sets, such as [`TD`](@ref).
```math
P_\\ell = \\frac{E_\\ell}{(\\sqrt{V_\\ell W_\\ell})^p}
```

**`accept_rate`** is the acceptance rate in the accept-reject method for chosing new indices from the active set for further refinement. A lower value means a more global search stretegy. Usefull when the adaptive method finds itself stuck along a single coordinate direction. Default value is 1 (no accept-reject).

# Keywords valid for [`QMC`](@ref) sample methods

**`nb_of_shifts`** indicates the number of randomly-shifted rank-1 lattice rules used in the [`QMC`](@ref) method. The higher this value, the more accurate the variance estimation will be. This keyword can also be a function that returns the number of shifts on each level (or index). Typically, the number of shifts on coarser levels will be smaller than the number of shifts on finer levels, since, on these coarse levels, many QMC samples will be available. Default is 10 random shifts.

!!! note

    We recommend a value for `nb_of_shifts` between 10 and 30.

**`sample_mul_factor`** is the multplication factor in the QMC "doubling algorithm". Typically, and as can be justified theoretically, this value is equal to 2. We find that it is sometimes more beneficial to use a smaller multiplication factor, since direclty doubling the number of samples is a quite dramatic event. Default: 1.2.

**`point_generator`** is a keyword that allows the user to specify a custom QMC point set of type `AbstractLatticeRule`. Default is a 3600-dimensional lattice rule from [this website](https://web.maths.unsw.edu.au/~fkuo/lattice/index.html). 

See also: [`run`](@ref), [`History`](@ref)
"""
struct Estimator{I<:AbstractIndexSet, S<:AbstractSampleMethod, T1, T2, T3}
    index_set::I
    sample_function::Function
    distributions::T1
    options::T2
    internals::T3
end

function Estimator(index_set::AbstractIndexSet, sample_method::AbstractSampleMethod, sample_function::Function, distributions::AbstractVector{<:AbstractDistribution}; kwargs...)

    # read options
    options = Dict{Symbol, Any}(kwargs)
    valid_options = get_valid_options(index_set, sample_method)
    for option in keys(options)
        if option ∉ valid_options
            throw(ArgumentError(string("in Estimator, invalid option ", option, " found")))
        end
    end
	options[:distributions] = distributions
    for option in valid_options
        parse!(index_set, sample_method, options, option)
    end

    # create estimator internals
    internals = EstimatorInternals(index_set, sample_method, options)

	Estimator{typeof(index_set), typeof(sample_method), typeof(distributions), typeof(options), typeof(internals)}(index_set, sample_function, distributions, options, internals)
end

"""
    Estimator(index_set::AbstractIndexSet,
              sample_method::AbstractSampleMethod,
              sample_function::Function,
              distribution::Distribution; kwargs...)

Same as [`Estimator`](@ref), but now with only one random parameter with specified distribution `distribution`.
"""
Estimator(index_set::AbstractIndexSet, sample_method::AbstractSampleMethod, sample_function::Function, distribution::AbstractDistribution; kwargs...) = Estimator(index_set, sample_method, sample_function, [distribution]; kwargs...)

for f in fieldnames(Estimator)
    @eval begin
        $f(estimator::Estimator) = estimator.$f
    end
end

## output formatting ##
show(io::IO, estimator::Estimator{I, S}) where {I, S} = print(io, string("Estimator{", index_set(estimator), ", ", S, "}"))
