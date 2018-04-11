# sample.jl : sample methods

## main sample function ##
function sample!(estimator::Estimator,index::Index,nb_of_samples::N where {N<:Integer}, warm_up::Bool)

    estimator.verbose && begin
        name = warm_up ? "warm-up" : "additional"
        println("Taking $(nb_of_samples) $(name) samples at level $(index[1])...")
    end

    # add new index if no samples are taken yet
    if warm_up
        estimator.nsamples[index] = 0
        estimator.samples[index] = eltype(values(estimator.samples))()
        estimator.total_work[index] = 0.
    end

    # parallel sampling
    istart = estimator.nsamples[index]+1
    iend = istart + nb_of_samples - 1
    parallel_sample!(estimator,index,istart,iend)
end

## print (optimal) number of samples ##
function print_number_of_samples(estimator,samples)
    println("Samples will be updated according to")
    n = 14
    nb = 2
    border = spaces(nb)
    level_name = ndims(estimator.method) == 1 ? "level" : "index"
    header = string(border,level_name,spaces(n-nb-length(level_name)))
    header = string(header,"N",spaces(n))
    println(repeat("-",29))
    println(header)
    println(repeat("-",29))
    for index in sort(collect(keys(samples)))
        index_str = "$(index)"
        str = string(border,index_str,spaces(n-nb-length(index_str)))
        str = string(str,@sprintf("%s",samples[index]))
        println(str)
    end
end
