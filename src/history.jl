## history.jl : a history object for bookkeeping
#
# A history object contains bookkeeping data for multilevel estimators.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2019

## History ##
mutable struct History{T}
    data::T
end

History() = History(Vector{Dict{Symbol, Any}}(undef, 0))

function push!(history::History, estimator::Estimator, tol::Real, elapsed::Real)
    h = Dict{Symbol, Any}()
    h[:type]          = typeof(estimator)
    h[:ndims]         = ndims(estimator)
	h[:name]          = estimator[:name]
	h[:folder]        = estimator[:folder]
    h[:elapsed]       = elapsed
    h[:tol]           = tol
    h[:index_set]     = keys(estimator)
    h[:mse]           = mse(estimator)
    h[:rmse]          = rmse(estimator)
    h[:mean]          = mean(estimator)
    h[:var]           = var(estimator)
    h[:varest]        = varest(estimator)
    h[:bias]          = bias(estimator)
    h[:E]             = apply(mean0, estimator)
    h[:V]             = apply(var0, estimator)
    h[:dE]            = apply(mean, estimator)
    h[:dV]            = apply(var, estimator)
    h[:W]             = apply(work, estimator)
    h[:W]             = apply(time, estimator)
    h[:α]             = α(estimator)
    h[:β]             = β(estimator)
    h[:γ]             = γ(estimator)
    h[:nb_of_samples] = copy(nb_of_samples(estimator)) 

    # save samples if required
	if estimator[:save_samples]
        h[:samples] = samples(estimator)
        h[:samples_diff] = samples_diff(estimator)
    end

    push!(history.data, h)

    # save history
    save(history)
end

apply(f::Function, estimator::Estimator) = Dict(i => f(estimator, i) for i in keys(estimator))

save(history::History) = @save joinpath(history[:folder], history[:name]) history

getindex(history::History, s::Symbol) = history.data[end][s]

getindex(history::History, i::Integer) = history.data[i]

length(history::History) = length(history.data)

show(io::IO, history::History) = print(io, "MultilevelEstimators.jl history file")
