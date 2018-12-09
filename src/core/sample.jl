# sample.jl : sample methods

## main sample! function ##
function sample!(estimator::Estimator, index::Index, n::Integer)

    warm_up = !contains_samples_at_index(estimator, index)
    print_sample!_header(estimator, index, n, warm_up)
    warm_up && add_index(estimator, index)

    Istart = nb_of_samples(estimator, index) + 1
    Iend = Istart + n - 1
    t = @elapsed parallel_sample!(estimator, index, Istart, Iend)
    update_nb_of_samples!(estimator, index, n)
    update_total_work!(estimator, index, n)
    update_total_time!(estimator, index, t)

    print_sample!_footer()
end

## sample for MC ##
function parallel_sample!(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index, Istart::Integer, Iend::Integer)

    m = nb_of_uncertainties(estimator, index)
    f(i) = estimator.sample_function(index, transform.(view(distributions(estimator), 1:m), rand(m)))
    all_workers = workers()
    worker_idcs = 1:min(nb_of_workers(estimator, index), nworkers())
    pool = CachingPool(all_workers[worker_idcs])
    batch_size = ceil(Int, (Iend-Istart+1)/length(worker_idcs))
    retry_delays = ExponentialBackOff(n=4)
    S = pmap(f, pool, Istart:Iend, batch_size=batch_size, retry_delays=retry_delays)
    clear!(pool)

    s_diff = first.(S)
    s = last.(S)

    append_samples!(estimator, index, s_diff, s, Iend-Istart+1)
end

function append_samples!(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index, s_diff, s, n::Integer)
    for n_qoi in 1:nb_of_qoi(estimator)
        append_samples_diff!(estimator, n_qoi, index, getindex.(s_diff, n_qoi))
        append_samples!(estimator, n_qoi, index, getindex.(s, n_qoi))
    end
end

function append_samples!(estimator::Estimator{<:MG, <:MC}, index::Index, s_diff, s, n::Integer)
    R = CartesianIndices(UnitRange.(zero(index), index))
    for n_qoi in 1:nb_of_qoi(estimator)
        acc = fill(zero(eltype(eltype(eltype(s)))), n)
        for idx in R
            append_samples_diff!(estimator, n_qoi, Index(idx), getindex.(getindex.(s_diff, Ref(idx+one(idx))), n_qoi))
            append_samples!(estimator, n_qoi, Index(idx), getindex.(getindex.(s, Ref(idx+one(idx))), n_qoi))
            acc += weight(estimator, Index(idx))*getindex.(getindex.(s_diff, Ref(idx+one(idx))), n_qoi)
        end	
        append_accumulator!(estimator, n_qoi, acc)
    end
end

## sample for QMC index sets ##
