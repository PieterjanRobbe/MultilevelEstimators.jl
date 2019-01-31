## sample.jl : collection of methods to add samples to an Estimator
#
# Methods to add samples to the Estimator.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for
# Multilevel Monte Carlo Methods (c) Pieterjan Robbe, 2019

function update_samples(estimator::Estimator, n_opt)
    for τ in keys(n_opt)
        n_due = n_opt[τ] - nb_of_samples(estimator, τ)
        n_due > 0 && sample!(estimator, τ, n_due)
    end
end

function update_samples(estimator::Estimator{<:U}, n_opt)
    for τ in keys(n_opt)
		n_opt[τ] > 0 && sample!(estimator, τ, n_opt[τ])
    end
end

function sample!(estimator::Estimator, index::Index, n::Integer)
    warm_up = !has_samples_at_index(estimator, index)
	estimator[:verbose] && print_sample!_header(estimator, index, n, warm_up)

    Istart = nb_of_samples(estimator, index) + 1
    Iend = Istart + n - 1
    t = @elapsed parallel_sample!(estimator, index, Istart, Iend)

    add_to_nb_of_samples(estimator, index, n)
    add_to_total_work(estimator, index, n)
    add_to_total_time(estimator, index, t)

	estimator[:verbose] && print_sample!_footer()
end

#
# MC
#
function parallel_sample!(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index, Istart::Integer, Iend::Integer)

    m = estimator[:nb_of_uncertainties](index)
    f(i) = estimator.sample_function(index, transform.(view(distributions(estimator), 1:m), rand(m)))
    all_workers = workers()
    worker_idcs = 1:min(estimator[:nb_of_workers](index), nworkers())
    pool = CachingPool(all_workers[worker_idcs])
    batch_size = ceil(Int, (Iend - Istart + 1)/length(worker_idcs))
    retry_delays = ExponentialBackOff(n = 3)
    S = pmap(f, pool, Istart:Iend, batch_size = batch_size, retry_delays = retry_delays)
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

function append_samples!(estimator::Estimator{<:U, <:MC}, index::Index, s_diff, s)
    for idx in zero(index):index
		for n_qoi in 1:estimator[:nb_of_qoi]
			append_samples_diff!(estimator, n_qoi, idx, getindex.(getindex.(s_diff, Ref(idx + one(idx))), n_qoi))
			append_samples!(estimator, n_qoi, idx, getindex.(getindex.(s, Ref(idx + one(idx))), n_qoi))
			idx == zero(idx) && append_to_accumulator!(estimator, n_qoi, zeros(length(s_diff)))
			add_to_accumulator!(estimator, n_qoi, 1 / Prob(estimator, idx) * getindex.(getindex.(s_diff, Ref(idx + one(idx))), n_qoi))
        end
    end
end

#
# QMC
#
function parallel_sample!(estimator::Estimator{<:AbstractIndexSet, <:QMC}, index::Index, Istart::Integer, Iend::Integer)

    m = estimator[:nb_of_uncertainties](index)
    x(i) = begin
		point = get_point(generator(estimator, index, i[1]), i[2] - 1)
		length(point) < m && append!(point, rand(m - length(point)))
		view(point, 1:m)
    end
    f(i) = estimator.sample_function(index, transform.(view(distributions(estimator), 1:m), x(i)))
    all_workers = workers()
    worker_idcs = 1:min(estimator[:nb_of_workers](index), nworkers())
    pool = CachingPool(all_workers[worker_idcs])
    batch_size = ceil(Int, (Iend - Istart + 1)/length(worker_idcs))
    retry_delays = ExponentialBackOff(n = 3)
    S = pmap(f, pool, Iterators.product(1:estimator[:nb_of_shifts](index), Istart:Iend), batch_size = batch_size, retry_delays = retry_delays)
    clear!(pool)

    s_diff = first.(S)
    s = last.(S)

    append_samples!(estimator, index, s_diff, s)
end

function append_samples!(estimator::Estimator{<:AbstractIndexSet, <:QMC}, index::Index, s_diff, s)
    for n_shift in 1:estimator[:nb_of_shifts](index)
        for n_qoi in 1:estimator[:nb_of_qoi]
            append_samples_diff!(estimator, n_qoi, n_shift, index, getindex.(view(s_diff, n_shift, :), n_qoi))
            append_samples!(estimator, n_qoi, n_shift, index, getindex.(view(s, n_shift, :), n_qoi))
        end
    end
end

function append_samples!(estimator::Estimator{<:U, <:QMC}, index::Index, s_diff, s)
    for idx in zero(index):index
        for n_shift in 1:estimator[:nb_of_shifts](index)
            for n_qoi in 1:estimator[:nb_of_qoi]
				append_samples_diff!(estimator, n_qoi, n_shift, idx, getindex.(getindex.(view(s_diff, n_shift, :), Ref(idx + one(idx))), n_qoi))
				append_samples!(estimator, n_qoi, n_shift, idx, getindex.(getindex.(view(s, n_shift, :), Ref(idx + one(idx))), n_qoi))
				idx == zero(idx) && append_to_accumulator!(estimator, n_qoi, n_shift, zeros(size(s_diff, 2)))
				add_to_accumulator!(estimator, n_qoi, n_shift, 1 / Prob(estimator, idx) * getindex.(getindex.(view(s_diff, n_shift, :), Ref(idx + one(idx))), n_qoi))
            end
        end
    end
end
