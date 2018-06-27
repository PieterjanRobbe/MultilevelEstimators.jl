## history.jl : store general info about an estimator

mutable struct History{N<:Integer}
    iter::N
    time_stamp::DateTime
    data::Vector{Dict{Symbol, Any}}
end

History() = History(0,now(),Dict{Symbol,Any}[])

function push!(h::History, estimator::Estimator,tol::T where {T<:Real})
    old_time_stamp = h.time_stamp

    h.iter += 1 # add extra iteration
    h.time_stamp = now() # update time stamp
    push!(h.data,Dict{Symbol,Any}()) # add empty Dict

    # log all keys
    h[:name] = estimator.name
    h[:tol] = tol
    h[:index_set] = keys(estimator)
    h[:folder] = estimator.folder
    h[:runtime] = Dates.value(h.time_stamp - old_time_stamp)/1e3 # milliseconds
    h[:cost] = sum([estimator.nb_of_shifts*estimator.nsamples[index]*cost(estimator,index) for index in keys(estimator)])
    h[:nsamples] = [estimator.nb_of_shifts*estimator.nsamples[index] for index in keys(estimator)]
    h[:nb_of_shifts] = estimator.nb_of_shifts
    h[:mse] = mse(estimator)
    h[:rmse] = rmse(estimator)
    h[:mean] = mean(estimator)
    h[:var] = var(estimator)
    h[:varest] = varest(estimator)
    h[:bias] = bias(estimator)
    h[:E] = [mean0(estimator,index) for index in keys(estimator)]
    h[:V] = [var0(estimator,index) for index in keys(estimator)]
    h[:dE] = [mean(estimator,index) for index in keys(estimator)]
    h[:dV] = [var(estimator,index) for index in keys(estimator)]
    h[:W] = [cost(estimator,index) for index in keys(estimator)]
    h[:α] = α(estimator)
    h[:β] = β(estimator)
    h[:γ] = γ(estimator)
    if estimator.store_samples
        h[:samples] = estimator.samples
    end
    h[:ndims] = ndims(estimator)
	h[:method] = estimator.method
	h[:adaptive_index_set] = estimator.adaptive_index_set

    # save
    save(h)
end

# save history
save(h::History) = 
jldopen(string(h[:folder],"history.jld"), "w") do file
    addrequire(file, MultilevelEstimators)
    write(file, "history", h)
end

# convenience functions
setindex!(h::History, val, s::Symbol) = h.data[h.iter][s] = val
getindex(h::History, s::Symbol) = h.data[h.iter][s]
getindex(h::History, i::N where {N<:Int}) = h.data[i]

# show methods
show(io::IO, history::History) = print(io, "MultilevelEstimators.jl history file")
