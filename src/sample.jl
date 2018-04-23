# sample.jl : sample methods

## main sample function ##
function sample!(estimator::Estimator,index::Index,nb_of_samples::N where {N<:Integer})

    warm_up = !haskey(estimator.samples[1],index)

    estimator.verbose && begin
        name = warm_up ? "warm-up" : "additional"
        is_qmc = isa(estimator,QuasiMonteCarloTypeEstimator)
        mult = is_qmc ? "$(estimator.nb_of_shifts) x " : ""
        withs = nb_of_samples > 1 ? "s" : ""
        println("Taking $(mult)$(nb_of_samples) $(name) sample$(withs) at level $(index[1])...")
    end

    # add new index if no samples are taken yet
    if warm_up
        estimator.nsamples[index] = 0
        for n_qoi = 1:estimator.nb_of_qoi
            for n_shift = 1:estimator.nb_of_shifts
                estimator.samples[n_qoi,n_shift][index] = valtype(estimator.samples[n_qoi])()
                estimator.samples0[n_qoi,n_shift][index] = valtype(estimator.samples0[n_qoi])()
            end
        end
        estimator.total_work[index] = 0.
        estimator.number_generators[index] = random_shift(estimator.number_generator)
    end

    # parallel sampling
    istart = estimator.nsamples[index] + 1
    iend = istart + nb_of_samples - 1
    estimator.parallel_sample_function(estimator,index,istart,iend)
end
