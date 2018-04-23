# run.jl : run estimator

function run(estimator::Estimator, tols::Vector{T} where {T<:AbstractFloat})

    # input checking
    all(tols.>0) || throw(ArgumentError("supplied tolerance(s) must be positive, got $(tols)"))

    # make history
    h = History()

    # run the sequence of tolerances
    for tol in tols
        _run(estimator,tol)
        push!(h,estimator,tol) # log the results in history
        clear(estimator) # prepare new run 
    end

    return h
end

function run(estimator::Estimator, tol::T where {T<:AbstractFloat})
    if estimator.continuate
        tols = estimator.p0.^(estimator.ntols-1:-1:0)*tol
    else
        tols = [tol]
    end
    run(estimator,tols)
end
