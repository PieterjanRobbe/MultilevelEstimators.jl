## run.jl : entry point method for simulating an Estimator
#
# Runs the Estimator for the given tolerance, and logs a History entry.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

"""
```julia
run(estimator, ϵ)
```
Run the estimator and compute the expected value of the quantity of interest up to the given tolerance(s) `ϵ`.

# Examples
```jldoctest
julia>

```
"""
function run(estimator::Estimator, tols::AbstractVector{<:Real})

    # input checking
	all(tols.>0) || throw(ArgumentError(string("supplied tolerance(s) must be positive, got ", tols)))

    # make history
    #h = History()

    # run the sequence of tolerances
    for tol in tols
        _run(estimator,tol)
        clear(estimator)
        #push!(h,estimator,tol) # log the results in history
        #clear(estimator) # prepare new run 
    end

    #return h
end

run(estimator::Estimator, tol::Real) = run(estimator, get_tols(estimator, tol))
