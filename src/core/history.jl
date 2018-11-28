## history.jl : a history object for bookkeeping
#
# A history object contains bookkeeping data for multilevel estimators.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## History ##
mutable struct History{T}
    t_start::DateTime
    data::T
end

History() = History(now(), Vector{Dict{Symbol, Any}}(undef, 0))

function push!(history::History, estimator::Estimator, tol::Real)
    h = Dict{Symbol, Any}()
    h[:type]          = typeof(estimator)
    h[:ndims]         = ndims(estimator.index_set)
    h[:name]          = name(estimator)
    h[:folder]        = folder(estimator)
    h[:time_stamp]    = now()
    h[:tol]           = tol
    h[:index_set]     = keys(estimator)
    h[:mse]           = mse(estimator)
    h[:rmse]          = rmse(estimator)
    h[:mean]          = mean(estimator)
    h[:var]           = var(estimator)
    h[:varest]        = varest(estimator)
    h[:bias]          = bias(estimator)
    h[:E]             = apply_f(estimator, mean0)
    h[:V]             = apply_f(estimator, var0)
    h[:dE]            = apply_f(estimator, mean)
    h[:dV]            = apply_f(estimator, var)
    h[:W]             = apply_f(estimator, cost)
    h[:α]             = α(estimator)
    h[:β]             = β(estimator)
    h[:γ]             = γ(estimator)
    h[:nb_of_samples] = copy(nb_of_samples(estimator)) 

    # save samples if required
    if save_samples(estimator)
        h[:samples] = samples(estimator)
        h[:samples_diff] = samples_diff(estimator)
    end

    push!(history.data, h)

    # save history
    save(history)
end

apply_f(estimator::Estimator, f::Function) = Dict(i => f(estimator, i) for i in keys(estimator))

## save history ##
save(history::History) = @save joinpath(history[:folder], history[:name]) history

## convenience functions ##
getindex(history::History, s::Symbol) = history.data[end][s]
getindex(history::History, i::Integer) = history.data[i]
length(history::History) = length(history.data)

## show ##
show(io::IO, history::History) = print(io, "MultilevelEstimators.jl history file")
