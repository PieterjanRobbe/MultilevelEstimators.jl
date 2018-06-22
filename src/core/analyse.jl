# analyse.jl : analyse estimator

function analyse(estimator::Estimator; nsamples=0::N where {N<:Int}, max_time=Inf::T where {T<:Real})

    # input checking
    nsamples >= 0 || throw(ArgumentError("nsamples must be greater than zero, got $(nsamples)"))
    nsamples = nsamples == 0 ? 1000 : nsamples # default is 1000
    max_time > 0 || throw(ArmentError("max_time must be positive, got $(max_time)"))

    # make history
    h = History()

    # repeat until enough samples are taken, or until max_time is exceeded
    done = false
    cntr = 0
    time = 0.
    while !done
        cntr += 1
        _analyse(estimator)
        push!(h,estimator,cntr) # log the results in history
        time = cumsum([h[i][:runtime] for i in 1:h.iter])[end]
        done = estimator.nsamples[first(keys(estimator))] >= nsamples || time > max_time
    end

    return h
end

function _analyse(estimator::Estimator)

    # first run: take 1 sample at each level/index
    if isempty(estimator.nsamples)
        for level in get_index_set(estimator.method,estimator.max_level)
            sample!(estimator,level,1)
            push!(estimator,level)
        end

    # all other runs: scale number of samples such that equal amount of work on each level/index
    else
        for level in keys(estimator)
            n = round(Int,1/cost(estimator,level)*cost(estimator,keys(estimator)[indmax([cost(estimator,level) for level in keys(estimator)])]))
            sample!(estimator,level,n)
        end
    end
    
    estimator.verbose && print_status(estimator)
end
