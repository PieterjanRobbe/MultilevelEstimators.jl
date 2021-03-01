## sample.jl : collection of methods to add samples to an Estimator
#
# Methods to add samples to the Estimator.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for
# Multilevel Monte Carlo Methods (c) Pieterjan Robbe, 2021

function update_samples(estimator::Estimator, n_opt)
    flag = false
    for τ in keys(n_opt)
        n_due = n_opt[τ] - nb_of_samples(estimator, τ)
        if n_due > 0
            ret = sample!(estimator, τ, n_due)
            if ret
                flag = true
            end
        end
    end
    return flag
end

function update_samples(estimator::Estimator{<:U}, n_opt)
    flag = false
    for τ in keys(n_opt)
        if n_opt[τ] > 0
            ret = sample!(estimator, τ, n_opt[τ])
            if ret
                flag = true
            end
        end
    end
    return flag
end

function sample!(estimator::Estimator{I}, index::Index, n::Integer) where I
    warm_up = !has_samples_at_index(estimator, index)
    estimator[:verbose] && print_sample!_header(estimator, index, n, warm_up)

    Istart = I <: U ? sum(values(nb_of_samples(estimator))) + 1 : nb_of_samples(estimator, index) + 1 # count zero level if U
    Iend = Istart + n - 1
    t, flag = parallel_sample(estimator, index, Istart, Iend)

    flag && return true

    add_to_nb_of_samples(estimator, index, n)
    add_to_total_work(estimator, index, n)
    add_to_total_time(estimator, index, t)

    estimator[:verbose] && print_sample!_footer()

    return false
end

function parallel_sample(estimator::Estimator, index::Index, Istart::Integer, Iend::Integer)
    if estimator[:checkpoint]
        t, flag, S = parallel_sample_checkpoint(estimator, index, Istart, Iend)
    else
        t = @elapsed flag, S = parallel_sample_online(estimator, index, Istart, Iend)
    end

    flag && return true

    s_diff = first.(S)
    s = last.(S)

    append_samples!(estimator, index, s_diff, s)

    return t, false
end

#
# MC
#
function parallel_sample_online(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index, Istart::Integer, Iend::Integer)

    m = estimator[:nb_of_uncertainties](index)
    f(i) = estimator.sample_function(index, transform.(view(distributions(estimator), 1:m), rand(m)))
    all_workers = workers()
    worker_idcs = 1:min(estimator[:nb_of_workers](index), nworkers())
    pool = CachingPool(all_workers[worker_idcs])
    #batch_size = ceil(Int, (Iend - Istart + 1)/length(worker_idcs))
    batch_size = 1
    retry_delays = ExponentialBackOff(n = 3)
    S = pmap(f, pool, Istart:Iend, batch_size = batch_size, retry_delays = retry_delays)
    clear!(pool)

    return false, S
end

function parallel_sample_checkpoint(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index, Istart::Integer, Iend::Integer)

    S = Vector{Tuple{Matrix{Float64}, Matrix{Float64}}}(undef, Iend - Istart + 1)
    restart = Dict(i => false for i in Istart:Iend)
	t = 0.0

    # loop over all samples
    for i in Istart:Iend
        dir = joinpath(estimator[:samples_dir], join(index.I, "_"), string(i))
        dQ_file = joinpath(dir, "dQ.dat")
        Qf_file = joinpath(dir, "Qf.dat")
		W_file = joinpath(dir, "time.dat")
        # check if sample already exists
        if isfile(dQ_file) && isfile(Qf_file)
            dQ = readdlm(dQ_file)
            Qf = readdlm(Qf_file)
            S[i - Istart + 1] = (dQ, Qf)
			if isfile(W_file) # read time of sample if time file exists
				t += readdlm(W_file)[1]
			elseif !(estimator[:cost_model] isa EmptyFunction) 
				str = string("File time.dat not found for sample number ", i, ", timings for ", print_elname(estimator), " ", index, " will be unreliable!") 
				@warn str
			end
        else # if not, write parameter values
            m = estimator[:nb_of_uncertainties](index)
            params = transform.(view(distributions(estimator), 1:m), rand(m))
            isdir(dir) || mkpath(dir)
            writedlm(joinpath(dir, "params.dat"), params)
            restart[i] = true
        end
    end

    estimator[:verbose] && print_restart(index, restart, estimator[:samples_dir], "*")

    # check for restart
    if any(values(restart))
        estimator[:verbose] && print_footer()
        restart_flag = true
    else
        restart_flag = false
    end

    return t, restart_flag, S
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
function parallel_sample_online(estimator::Estimator{I, <:QMC}, index::Index, Istart::Integer, Iend::Integer) where I<:AbstractIndexSet

    m = estimator[:nb_of_uncertainties](index)
    x(i) = begin
        _index = I <: U ? zero(index) : index # take points from level 0 if unbiased
        point = get_point(generator(estimator, index, i[1]), i[2] - 1)
        length(point) < m && append!(point, rand(m - length(point)))
        view(point, 1:m)
    end
    f(i) = estimator.sample_function(index, transform.(view(distributions(estimator), 1:m), x(i)))
    all_workers = workers()
    worker_idcs = 1:min(estimator[:nb_of_workers](index), nworkers())
    pool = CachingPool(all_workers[worker_idcs])
    #batch_size = ceil(Int, (Iend - Istart + 1)/length(worker_idcs))
    batch_size = 1
    retry_delays = ExponentialBackOff(n = 3)
    S = pmap(f, pool, Iterators.product(1:estimator[:nb_of_shifts](index), Istart:Iend), batch_size = batch_size, retry_delays = retry_delays)
    clear!(pool)

    return false, S
end

function parallel_sample_checkpoint(estimator::Estimator{I, <:QMC}, index::Index, Istart::Integer, Iend::Integer) where I<:AbstractIndexSet

    K = estimator[:nb_of_shifts](index)
    S = Matrix{Tuple{Matrix{Float64}, Matrix{Float64}}}(undef, K, Iend - Istart + 1)
    restart = Dict(i => false for i in Istart:Iend)

    # loop over all samples
    for k in 1:K
        for i in Istart:Iend
            dir = joinpath(estimator[:samples_dir], join(index.I, "_"), string(i), string(k))
            dQ_file = joinpath(dir, "dQ.dat")
            Qf_file = joinpath(dir, "Qf.dat")
            # check if sample already exists
            if isfile(dQ_file) && isfile(Qf_file)
                dQ = readdlm(dQ_file)
                Qf = readdlm(Qf_file)
                S[k, i] = (dQ, Qf)
            else # if not, write parameter values
                m = estimator[:nb_of_uncertainties](index)
                _index = I <: U ? zero(index) : index
                point = get_point(generator(estimator, index, k), i - 1)
                length(point) < m && append!(point, rand(m - length(point)))
                params = transform.(view(distributions(estimator), 1:m), view(point, 1:m))
                isdir(dir) || mkpath(dir)
                writedlm(joinpath(dir, "params.dat"), params)
                restart[i] = true
            end
        end
    end

    estimator[:verbose] && print_restart(index, restart, estimator[:samples_dir], joinpath("*", "*"))

    # check for restart
    if any(values(restart))
        estimator[:verbose] && print_footer()
        restart_flag = true
    else
        restart_flag = false
    end

    return restart_flag, S
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
