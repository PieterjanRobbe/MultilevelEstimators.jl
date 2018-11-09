## multigrid_multilevel_quasi_monte_carlo.jl : run Multilevel Quasi-Monte Carlo estimator

## MultiGrid Multilevel Quasi-Monte Carlo parallel sampling ##
function parallel_sample!(estimator::MultiGridMultiLevelQuasiMonteCarloEstimator,index::Level,istart::N,iend::N) where {N<:Integer}

    # parallel sampling
    wp = CachingPool(workers())
    f(i) = estimator.sample_function(index,get_point(estimator.number_generators[index],i[2],i[1]),estimator.user_data)
    nshifts = estimator.nb_of_shifts
    nqoi = estimator.nb_of_qoi
    t = @elapsed all_samples = pmap(wp,f,Base.Iterators.product(1:nshifts,istart:iend),batch_size=ceil(Int,(iend-istart+1)/nworkers()), retry_delays = ExponentialBackOff(n = 3))

    # extract samples
    samples = last.(all_samples)
    dsamples = first.(all_samples)

    # append samples
    for idx in Iterators.product(UnitRange.(0,index.+1)...)
        for j in 1:nshifts
            for i in 1:nqoi
                append!(estimator.samples[i,j][idx],getindex.(getindex.(dsamples[j,:],(idx.+1)...),i))
                append!(estimator.samples0[i,j][idx],getindex.(getindex.(samples[j,:],(idx.+1)...)i))
            end
        end
    end
    estimator.nsamples[index] += iend-istart+1
    estimator.total_work[index] += estimator.use_cost_model ? (iend-istart+1)*estimator.cost_model(index) : t
end
