# sample.jl : sample methods

## main sample function ##
function sample!(estimator::Estimator,index::Index,nb_of_samples::N where {N<:Integer})

    warm_up = !haskey(estimator.samples,index)

    estimator.verbose && begin
        name = warm_up ? "warm-up" : "additional"
        println("Taking $(nb_of_samples) $(name) samples at level $(index[1])...")
    end

    # add new index if no samples are taken yet
    if warm_up
        estimator.nsamples[index] = 0
        estimator.samples[index] = eltype(values(estimator.samples))()
        estimator.samples0[index] = eltype(values(estimator.samples0))()
        estimator.total_work[index] = 0.
    end

    # parallel sampling
    istart = estimator.nsamples[index]+1
    iend = istart + nb_of_samples - 1
    estimator.parallel_sample_function(estimator,index,istart,iend)
end
