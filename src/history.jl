## history.jl : a history object for bookkeeping
#
# A history object contains bookkeeping data for multilevel estimators.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for
# Multilevel Monte Carlo Methods (c) Pieterjan Robbe, 2019

## History ##
"""
    History

A History is the result of a call to [`run`](@ref) and contains useful diagnostic information about the state of the `Estimator`.

The information is stored as a `Vector` of `Dict`'s, where each entry in the vector is the result of a single run. See below for a list of keys added to these `Dict`s. Access to the key corresponding to the last entry in the `History` (that is, the one with the highest accuracy) is provided as `history[key]`.

See also the package [`Reporter`](https://github.com/PieterjanRobbe/Reporter.jl) for automatic generation of reports based on a History object.

# Keys

- `type::AbstractIndexSet`: the index set type of the `Estimator`.
- `ndims::Integer`: the number of dimensions of the `Estimator`.
- `name::String`: the name of this `Estimator`.
- `folder::String`: the directory where the results are stored.
- `elapsed::Real`: the elapsed time (in seconds).
- `tol::Real`: the request accuracy.
- `current_index_set::Vector{::Index}`: the levels (or indices) currently in use by the `Estimator`.
- `index_set::Vector{::Index}`: all levels (or indices) that are added to the `Estimator`.
- `mse::Real`: the measured mean square error of the `Estimator`. 
- `rmse::Real`: the measured root mean square error of the `Estimator`. 
- `mean::Real`: the estimated expected value of the quantity of interest, up to an absolute error of `tol`.
- `var::Real`: the measured variance of the quantity of interest.
- `varest::Real`: the measured variance of the `Estimator` for the expected value of the quantity of interest.
- `bias::Real`: the measured bias of the `Estimator`.
- `E::Dict{::Index, ::Real}`: a `Dict` with the measured expected value of the quantity of interest on each level (or index).
- `V::Dict{::Index, ::Real}`: a `Dict` with the measured variance of the quantity of interest on each level (or index).
- `dE::Dict{::Index, ::Real}`: a `Dict` with the measured expected value of the **difference** of the quantity of interest on each level (or index).
- `dV::Dict{::Index, ::Real}`: a `Dict` with the measured variance of the **difference** of the quantity of interest on each level (or index).
- `T::Dict{::Index, ::Real}`: a `Dict` with the measured computation time for a sample of the quantity of interest on each level (or index).
- `W::Dict{::Index, ::Real}`: a `Dict` with the computational cost for a sample of the quantity of interest on each level (or index). *Only available when a cost model was provided!*
- `α::Tuple`: the estimated rates of decay of the expected value of the difference of the quantity of interest.
- `β::Tuple`: the estimated rates of decay of the variance of the difference of the quantity of interest.
- `γ::Tuple`: the estimated rates of increase of the cost of the difference of the quantity of interest.
- `nb_of_samples::Array`: an `Array` with, on each index, the number of samples taken on that index. 
- `logbook`: a logbook with information about the adaptive algorithm *Only available for `Estimator`s with index sets of type `AD`*
- `samples::Array`: an `Array` with on each index, the collection of samples of the quantity of interest taken on that index. *Only available when the optional key `save_samples = true` was used* 
- `samples_diff::Array`: an `Array` with on each index, the collection of samples of the **difference** of the quantity of interest taken on that index. *Only available when the optional key `save_samples = true` was used* 
See also: [`Estimator`](@ref), [`run`](@ref)
"""
mutable struct History{T}
    data::T
end

History() = History(Vector{Dict{Symbol, Any}}(undef, 0))

function push!(history::History, estimator::Estimator, tol::Real, elapsed::Real)
    h = Dict{Symbol, Any}()
    h[:type]              = typeof(estimator)
    h[:ndims]             = ndims(estimator)
	h[:name]              = estimator[:name]
	h[:folder]            = estimator[:folder]
    h[:elapsed]           = elapsed
    h[:tol]               = tol
	h[:current_index_set] = copy(collect(keys(estimator)))
	h[:index_set]         = copy(collect(all_keys(estimator)))
    h[:mse]               = mse(estimator)
    h[:rmse]              = rmse(estimator)
    h[:mean]              = mean(estimator)
    h[:var]               = var(estimator)
    h[:varest]            = varest(estimator)
    h[:bias]              = bias(estimator)
    h[:E]                 = apply(mean0, estimator)
    h[:V]                 = apply(var0, estimator)
    h[:dE]                = apply(mean, estimator)
    h[:dV]                = apply(var, estimator)
    h[:T]                 = apply(time, estimator)
    h[:α]                 = α(estimator)
    h[:β]                 = β(estimator)
    h[:γ]                 = γ(estimator)
    h[:nb_of_samples]     = copy(nb_of_samples(estimator)) 

	# add cost model is it was provided
	if !(estimator[:cost_model] isa EmptyFunction)
    	h[:W]             = apply(cost, estimator)
	end

	# add logbook when using adpative index sets
	if estimator isa Estimator{<:AD}
		h[:logbook]       = copy(logbook(estimator))
	end

    # save samples if required
    if estimator[:save]
        if estimator[:save_samples]
            h[:samples] = samples(estimator)
            h[:samples_diff] = samples_diff(estimator)
        end

        push!(history.data, h)

        # save history
        save(history)
    end
end

apply(f::Function, estimator::Estimator) = Dict(i => f(estimator, i) for i in all_keys(estimator))

save(history::History) = @save joinpath(history[:folder], history[:name]) history

getindex(history::History, s::Symbol) = history.data[end][s]

getindex(history::History, i::Integer) = history.data[i]

haskey(history::History, s::Symbol) = haskey(history.data[end], s)

length(history::History) = length(history.data)

show(io::IO, history::History) = print(io, "MultilevelEstimators.jl history file")
