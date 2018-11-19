# sample.jl : sample methods

## main sample! function ##
function sample!(estimator::Estimator, index::Index, nb_of_samples::Integer)

    warm_up = contains_samples_at_index(estimator, index)
    print_sample!_header(estimator, index, nb_of_samples, warm_up)
    warm_up && add_index(estimator, index)

    Istart = nb_of_samples_at_index(estimator, index) + 1
    Iend = Istart + nb_of_samples - 1
    parallel_sample!(estimator, index, Istart, Iend)
    update_nb_of_samples!(estimator, nb_of_samples)
    update_total_work!(estimator, index, t)
end

function parallel_sample(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index, Istart::Integer, Iend::Integer)

    f(i) = estimator.sample_function(index, transform.(distributions(estimator), rand(stochastic_dim(estimator))))
    wp = CachingPool(workers())
    batch_size = ceil(Int, (Iend-Istart+1)/nb_of_workers_at_index(estimator, index))
    retry_delays = ExponentialBackOff(n = 4)
    t = @elapsed S = pmap(f, wp, batch_size=batch_size, retry_delays=retry_delays)

    samples_diff = first.(D)
    samples = last.(S)

    for n_qoi in 1:nb_of_qoi(estimator)
        append!(samples_diff(estimator, n_qoi, index), samples_diff[n_qoi])
        append!(samples(estimator, n_qoi, index), samples[n_qoi])
    end
end
