## parse_sampler.jl : parse Sampler from input Dict

# setup will parse all inputs from the Dict and set the appropriate sampler settings
function sampler_setup{S<:AbstractString}(settings::Dict{S,Any})

    ## required ##

    # index_set
    if !haskey(settings,"index_set")
        error("I need an index_set (SL/ML/FT/TD/HC/AD) before I can do anything!")
    else
        index_set = settings["index_set"]
        delete!(settings,"index_set")
        if !(typeof(index_set) <: IndexSet)
            error("incorrect index_set specified!")
        end
    end
    d = ndims(indexSet)

    # number_generator
    if !haskey(settings,"number_generator")
        error("I need a number_generator before I can do anything!")
    else
        number_generator = settings["number_generator"]
        delete!(settings,"number_generator")
        if typeof(number_generator <: NumberGenerator)
            number_generators = Dict{Index{d,Vector{Int64}},typeof(number_generator)}()
            number_generators[zero_index(d)] = number_generator
        elseif ( typeof(number_generator) <: Dict && Base.valtype(number_generator) <: NumberGenerator )
            number_generators = number_generator
            if isempty(number_generators) || !haskey(number_generators,zero_index(d))
                error("number generator is empty or does not contain the zero $(d == 1 ? "level" : "index")")
            end
            generator_type = typeof(number_generators[zero_index(d)])
            for generator in number_generators
                if typeof(generator) != generator_type
                    error("number generators must be all of the same type!")
                end
            end
        else
            error("incorrect number_generator specified!")
        end
    end

    # sample_function
    if !haskey(settings,"sample_function")
        error("I need a sample_function (or QOI) before I can do anything!")
    else
        sampleFunction = settings["sample_function"]
        delete!(settings,"sample_function")
        if !(typeof(sample_function) <: Function)
            error("incorrect sample_function specified!")
        end
    end

    ## extra options ##

    if haskey(settings,"isMultigrid")
        isMultigrid = settings["isMultigrid"]
        if !(typeof(isMultigrid) <: Bool)
            throw(ArgumentError("isMultigrid must be of type Bool"))
        end
        delete!(settings,"isMultigrid")
        if haskey(settings,"isPrashant")
            isPrashant = settings["isPrashant"]
            if !(typeof(isPrashant) <: Bool)
                throw(ArgumentError("isMultigrid must be of type Bool"))
            end
            delete!(settings,"isPrashant")
        else
            isPrashant = false
        end
        if haskey(settings,"isCauchySchwarz")
            is_cauchy_schwarz = settings["isCauchySchwarz"]
            if !(typeof(is_cauchy_schwarz) <: Bool)
                throw(ArgumentError("isCauchySchwarz must be of type Bool"))
            end
            delete!(settings,"isCauchySchwarz")
        else
            is_cauchy_schwarz = false
        end
        if haskey(settings,"useBatches")
            use_batches = settings["useBatches"]
            if !(typeof(use_batches) <: Bool)
                throw(ArgumentError("useBatches must be of type Bool"))
            end
            delete!(settings,"useBatches")
        else
            use_batches = false
        end
        if haskey(settings,"giles_multigrid")
            giles_multigrid = settings["giles_multigrid"]
            if !(typeof(giles_multigrid) <: Bool)
                throw(ArgumentError("giles_multigrid must be of type Bool"))
            end
            delete!(settings,"giles_multigrid")
        else
            giles_multigrid = false
        end
    else
        isMultigrid = false
        isPrashant = false
        is_cauchy_schwarz = false
        use_batches = false
        giles_multigrid = false
    end
    ml_sample_fun = isMultigrid ? ( isPrashant ? prashant_sample_mg : ( giles_multigrid ? sample_gmg : sample_mg ) ) : sample

    if haskey(settings,"reuseSamples")
        reuseSamples = settings["reuseSamples"]
        if !(typeof(reuseSamples) <: Bool)
            throw(ArgumentError("reuseSamples must be of type Bool"))
        end
        delete!(settings,"reuseSamples")
    else
        reuseSamples = false
    end
    if reuseSamples && !isMultigrid
        throw(ArgumentError("for sample re-use, toggle multigrid on"))
    end

    if haskey(settings,"maxL")
        maxL = settings["maxL"]
        delete!(settings,"maxL")
        if typeof(maxL) <: Integer && maxL > 0
            mmaxL = maxL*ones(Int64,d)
        elseif typeof(maxL) == Vector{Int64} && all(maxL.>0)
            mmaxL = maxL
            maxL = mmaxL[1]
        else
            error("incorrect maxL specified!")
        end
    else
        maxL = 6
        mmaxL = maxL*ones(Int64,d)
    end
    if ( typeof(indexSet) == SL && maxL != 0 )
        warn("for SL simulation, maxL must equal zero. I will try to proceed with maxL=0")
        mmaxL = [0]
    end

    if haskey(settings,"costModel")
        icostModel = settings["costModel"]
        if !(typeof(icostModel) <: Function)
            error("incorrect cost model costModel specified!")
        end
        costModel = (index) -> mi_cost(icostModel,index)
    else
        costModel = (index) -> prod(2.^(index.indices))
    end

    if haskey(settings,"Z")
        Z = settings["Z"]
        delete!(settings,"Z")
        if !(typeof(Z) == Int64)
            error("incorrect number of qoi's Z specified!")
        elseif Z < 1
            error("number of qoi's Z must be positive!")
        end
    else
        Z = 1
    end

    if haskey(settings,"Nstar")
        Nstar = settings["Nstar"]
        delete!(settings,"Nstar")
        if !(typeof(Nstar) == Int64)
            error("incorrect number of warm-up samples Nstar specified!")
        elseif Nstar < 1
            error("number of warm-up samples Nstar must be positive!")
        end
    else
        Nstar = convert(Int64,ceil(32/nshifts(numberGenerators[Index(zeros(Int64,d))])))
    end

    if haskey(settings,"gaussianFieldSampler")
        gaussianFieldSampler = settings["gaussianFieldSampler"]
        delete!(settings,"gaussianFieldSampler")
        if !(typeof(gaussianFieldSampler) <: GaussianFieldSampler || (typeof(gaussianFieldSampler) <: Array && eltype(gaussianFieldSampler) <: GaussianFieldSampler ) )
            error("incorrect Gaussian field sampler gaussianFieldSampler specified!")
            # TODO need a better way to specify maxL for KL expansion and sampler, wrt adaptivity...
        elseif typeof(gaussianFieldSampler) == KLExpansion && maxL+1 < length(gaussianFieldSampler.eigenfunc)
            warn("you asked for $(maxL+1) levels, but the KL expansion only provides $(length(gaussianFieldSampler.eigenfunc)) eigenfunctions. I will try to continue anyway.")
        elseif typeof(gaussianFieldSampler) <: Array
            for i in 1:length(gaussianFieldSampler)
                if typeof(gaussianFieldSampler[i]) == KLExpansion && maxL+1 < length(gaussianFieldSampler[i].eigenfunc)
                    warn("you asked for $(maxL+1) levels, but the $(i)-th KL expansion only provides $(length(gaussianFieldSampler.eigenfunc)) eigenfunctions. I will try to continue anyway.")
                end
            end
        end
    else
        gaussianFieldSampler = EmptySampler()
    end

    if haskey(settings,"useTime")
        useTime = settings["useTime"]
        delete!(settings,"useTime")
        if !(typeof(useTime) == Bool)
            error("incorrect boolean useTime specified!")
        end
    else
        useTime = false
    end
    # check that either useTime is false and γ was provided or useTime is true and a cost model γ is provided
    if useTime == false && !haskey(settings,"costModel")
        warn("useTime is set to false, but no cost model costModel was provided. I will use the standard cost model and try to continue.")
    elseif useTime == true && haskey(settings,"costModel")
        warn("useTime is set to true, but also a cost model costModel was provided. I will use the true simulation times and try to continue.")
    end
    delete!(settings,"costModel")

    if haskey(settings,"safety")
        safety = settings["safety"]
        delete!(settings,"safety")
        if !(typeof(safety) == Bool)
            error("incorrect boolean safety specified!")
        end
    else
        safety = true
    end
    if Base.valtype(numberGenerators) <: QMCgenerator && safety == false
        warn("when using a QMC generator, safety must be turned on!")
        safety = true
    end

    if haskey(settings,"continuate")
        continuate = settings["continuate"]
        delete!(settings,"continuate")
        if !(typeof(continuate) == Bool)
            error("incorrect boolean continuate specified!")
        end
    else
        continuate = false
    end

    if haskey(settings,"nTOL")
        nTOL = settings["nTOL"]
        delete!(settings,"nTOL")
        if !(typeof(nTOL) == Int64)
            error("incorrect number of continuation tolerances nTOL specified!")
        elseif nTOL < 0
            error("number of continuation tolerances nTOL must be positive!")
        elseif !continuate
            warn("number of continuation tolerances nTOL was provided, but continuation is set to false. I will ignore entry 'nTOL'")
        end
    else
        nTOL = 10
    end

    if haskey(settings,"k")
        k = settings["k"]
        delete!(settings,"k")
        if !(typeof(k) == Tuple{Float64,Float64})
            error("incorrect trust parameters k specified!")
        elseif length(k) != 2
            error("I can only accept 2 trust parameters k0 and k1")
        elseif k[1]<0 || k[2]<0
            error("trust parameters k must be positive")
        elseif !continuate
            warn("trust parameters k where provided, but continuation is set to false. I will ignore entry 'k'")
        end
    else
        k = (0.1,0.1)
    end

    if haskey(settings,"showInfo")
        showInfo = settings["showInfo"]
        delete!(settings,"showInfo")
        if !(typeof(showInfo) == Bool)
            error("incorrect boolean showInfo specified!")
        end
    else
        showInfo = true
    end
    @debug showInfo = true

    if haskey(settings,"ioStream")
        ioStream = settings["ioStream"]
        if !(typeof(ioStream) == IO)
            error("incorrect ioStream specified!")
        end
    else
        ioStream = STDOUT
    end

    if haskey(settings,"storeSamples0")
        storeSamples0 = settings["storeSamples0"]
        delete!(settings,"storeSamples0")
        if !(typeof(storeSamples0) == Bool)
            error("incorrect boolean storeSamples0 specified!")
        end
    else
        storeSamples0 = false
    end

    if haskey(settings,"procMap")
        procMap = settings["procMap"]
        delete!(settings,"procMap")
        if !(typeof(procMap) == Dict{Index{d,Vector{Int64}},Int64})
            error("incorrect procMap specified!")
        elseif minimum(values(procMap)) <= 0
            error("wrong procMap provided, nprocs must be positive!")
        elseif maximum(values(procMap)) > nprocs()
            error("I only have $(nprocs()) processor(s) available. Lower the requested number of processors in the procMap or start Julia with more processors.")
        end
    else
        procMap = Dict{Index{d,Vector{Int64}},Int64}()
        dummyIndexSet = typeof(indexSet) <: AD ? FT(d) : indexSet
        id_x = getIndexSet(dummyIndexSet,max(20,maxL)) # hack
        [setindex!(procMap,nprocs(),i) for i in id_x]
    end

    if haskey(settings,"userType")
        userType = settings["userType"]
        delete!(settings,"userType")
    else
        userType = Void()
    end

    if haskey(settings,"max_indexset")
        max_indexset = settings["max_indexset"]
        delete!(settings,"max_indexset")
        if !(typeof(max_indexset) == Set{Index{d,Vector{Int64}}}) || isempty(max_indexset)
            error("Incorrect max_indexset specified!")
        end
    else
        inner = Set{Index{d,Vector{Int64}}}()
        for i in Base.product([1:mmaxL[i]-1 for i = 1:d]...)
            push!(inner,Index(i...))
        end
        outer = Set{Index{d,Vector{Int64}}}()
        for i in Base.product([1:mmaxL[i] for i = 1:d]...)
            push!(inner,Index(i...))
        end
        max_indexset = setdiff(outer, inner)
        delete!(max_indexset,mmaxL)
        push!(max_indexset,Index(mmaxL-1))
    end

    if !isempty(settings)
        str = ""
        for s in keys(settings)
            str *= s*", "
        end
        str = str[1:end-2]
        error("could not parse the following entries: "*str)
    end

    generatorStates = Dict{Index{d,Vector{Int64}},Int64}()

    samples = Dict{Index{d,Vector{Int64}},Array{Float64,3}}()
    samples0 = Dict{Index{d,Vector{Int64}},Array{Float64,3}}()

    T = Dict{Index{d,Vector{Int64}},Float64}()
    E = Dict{Index{d,Vector{Int64}},Vector{Float64}}()
    V = Dict{Index{d,Vector{Int64}},Vector{Float64}}()
    Vf = Dict{Index{d,Vector{Int64}},Vector{Float64}}()
    Wst = Dict{Index{d,Vector{Int64}},Float64}()
    W = Dict{Index{d,Vector{Int64}},Float64}()
    Vest = Dict{Index{d,Vector{Int64}},Vector{Float64}}()
    P = Dict{Index{d,Vector{Int64}},Float64}()

    return Sampler{d,typeof(indexSet),typeof(numberGenerators),typeof(gaussianFieldSampler),typeof(generatorStates),typeof(samples),typeof(T),typeof(E),typeof(userType),typeof(max_indexset)}(
                                                                                                                                                                                               indexSet, numberGenerators, sampleFunction, mmaxL, costModel, Z, Nstar, gaussianFieldSampler, useTime,
                                                                                                                                                                                               safety, continuate, nTOL, k, showInfo, ioStream, storeSamples0, procMap, userType, max_indexset, 
                                                                                                                                                                                               generatorStates, samples, samples0, T,E,V,Vf,Wst,W,Vest,P, ml_sample_fun,reuseSamples, Dict{Index{d,Vector{Int64}},Int64}(),is_cauchy_schwarz, use_batches, giles_multigrid)
end


