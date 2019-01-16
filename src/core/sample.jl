# sample.jl : sample methods

function sample!(estimator::Estimator, index::Index, n::Integer)
    warm_up = !has_samples_at_index(estimator, index)
    print_sample!_header(estimator, index, n, warm_up)

    Istart = nb_of_samples(estimator, index) + 1
    Iend = Istart + n - 1
    t = @elapsed parallel_sample!(estimator, index, Istart, Iend)

    add_to_nb_of_samples(estimator, index, n)
    add_to_total_work(estimator, index, n)
    add_to_total_time(estimator, index, t)
end

function parallel_sample!(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index, Istart::Integer, Iend::Integer)

	m = estimator[:nb_of_uncertainties](index)
    f(i) = estimator.sample_function(index, transform.(view(distributions(estimator), 1:m), rand(m)))
    all_workers = workers()
	worker_idcs = 1:min(estimator[:nb_of_workers](index), nworkers())
    pool = CachingPool(all_workers[worker_idcs])
    batch_size = ceil(Int, (Iend-Istart+1)/length(worker_idcs))
    retry_delays = ExponentialBackOff(n=4)
    S = pmap(f, pool, Istart:Iend, batch_size=batch_size, retry_delays=retry_delays)
    clear!(pool)

    s_diff = first.(S)
    s = last.(S)

    append_samples!(estimator, index, s_diff, s)
end

function append_samples!(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index, s_diff, s)
	for n_qoi in 1:estimator[:nb_of_qoi]
        append_samples_diff!(estimator, n_qoi, index, getindex.(s_diff, n_qoi))
        append_samples!(estimator, n_qoi, index, getindex.(s, n_qoi))
    end
end

function append_samples!(estimator::Estimator{<:MG, <:MC}, index::Index, s_diff, s)
    for idx in CartesianIndices(UnitRange.(zero(index), index))
        for n_qoi in 1:nb_of_qoi(estimator)
            append_samples_diff!(estimator, n_qoi, idx, getindex.(get_index.(s_diff, idx+one(idx)), n_qoi))
            append_samples!(estimator, n_qoi, idx, getindex.(getindex.(s, idx+one(idx)), n_qoi))
        end
    end
end
