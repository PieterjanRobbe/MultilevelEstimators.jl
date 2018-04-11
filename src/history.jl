## history.jl : store general info about an estimator

mutable struct History{N<:Integer}
    iter::N
    time_stamp::DateTime
    data::Vector{Dict{Symbol, Any}}
end

History() = History(0,now(),Dict{Symbol,Any}[])

function push!(h::History, estimator::Estimator)
    old_time_stamp = h.time_stamp

    h.iter += 1 # add extra iteration
    h.time_stamp = now() # update time stamp
    push!(h.data,Dict{Symbol,Any}()) # add empty Dict

    # log all keys
    h[:folder] = estimator.folder
    h[:runtime] = h.time_stamp - old_time_stamp 
    h[:nsamples] = estimator.nsamples
    h[:mse] = mse(estimator)
    h[:rmse] = rmse(estimator)
    h[:mean] = mean(estimator)
    h[:var] = var(estimator)
    if estimator.store_samples
        h[:samples] = estimator.samples
    end

    # save
    save(h)

    # clear estimator
    clear(estimator) # prepare new run 
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
