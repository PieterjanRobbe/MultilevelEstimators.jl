## estimator.jl : main Estimator type
#
# Representation of the main Estimator type.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2019

## Estimator ##
"""
    Estimator(index_set, sample_method, sample_function, distributions; kwargs...)

Create an `Estimator` with index set `index_set` and sample method `sample_method` for the expected value of the quantity of interest returned by `sample_function`, and where `distributions` is the uncertainty on the input parameters.

# Examples
```jldoctest
julia>

```
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

Estimator(index_set::AbstractIndexSet, sample_method::AbstractSampleMethod, sample_function::Function, distribution::AbstractDistribution; kwargs...) = Estimator(index_set, sample_method, sample_function, [distribution]; kwargs...)

for f in fieldnames(Estimator)
    @eval begin
        $f(estimator::Estimator) = estimator.$f
    end
end

ndims(::Estimator{<:AbstractIndexSet{d}}) where d = d 

get_index_set(estimator::Estimator, sz) = get_index_set(estimator.index_set, sz)

get_tols(estimator::Estimator, tol::T) where T<:Real = continuate(estimator) ? continuation_mul_factor(estimator).^(nb_of_tols(estimator)-1:-1:0)*tol : T[tol] 

## output formatting ##
show(io::IO, estimator::Estimator{I, S}) where {I, S} = print(io, string("Estimator{", index_set(estimator), ", ", S, "}"))
