## SPDE_init.jl : init functions for lognormal diffusion problem

## export statements ##
export init_SPDE_analyse_ml
export init_SPDE_analyse_mi
export init_SPDE_mc
export init_SPDE_mc_multiple
export init_SPDE_mlmc
export init_SPDE_mlmc_multiple
export init_SPDE_qmc
export init_SPDE_qmc_multiple
export init_SPDE_mlqmc
export init_SPDE_mlqmc_multiple
export init_SPDE_mimc
export init_SPDE_mimc_multiple
export init_SPDE_miqmc
export init_SPDE_miqmc_multiple
export init_SPDE_amimc
export init_SPDE_amimc_multiple
export init_SPDE_amiqmc
export init_SPDE_amiqmc_multiple
export init_SPDE_mgmlmc
export init_SPDE_mgmlmc_multiple

## init functions ##
init_SPDE_analyse_ml(;kwargs...) = init_SPDE(ML(),false,false,true,false;kwargs...)
init_SPDE_analyse_mi(;kwargs...) = init_SPDE(TD(2),false,false,true,false;kwargs...)
init_SPDE_mc(;kwargs...) = init_SPDE(SL(),false,false,false,false;kwargs...)
init_SPDE_mc_multiple(;kwargs...) = init_SPDE(SL(),false,true,false,false;kwargs...)
init_SPDE_mlmc(;kwargs...) = init_SPDE(ML(),false,false,false,false;kwargs...)
init_SPDE_mlmc_multiple(;kwargs...) = init_SPDE(ML(),false,true,false,false;kwargs...)
init_SPDE_qmc(;kwargs...) = init_SPDE(SL(),true,false,false,false;kwargs...)
init_SPDE_qmc_multiple(;kwargs...) = init_SPDE(SL(),true,true,false,false;kwargs...)
init_SPDE_mlqmc(;kwargs...) = init_SPDE(ML(),true,false,false,false;kwargs...)
init_SPDE_mlqmc_multiple(;kwargs...) = init_SPDE(ML(),true,true,false,false;kwargs...)
init_SPDE_mimc(;kwargs...) = init_SPDE(TD(2),false,false,false,false;kwargs...)
init_SPDE_mimc_multiple(;kwargs...) = init_SPDE(TD(2),false,true,false,false;kwargs...)
init_SPDE_miqmc(;kwargs...) = init_SPDE(TD(2),true,false,false,false;kwargs...)
init_SPDE_miqmc_multiple(;kwargs...) = init_SPDE(TD(2),true,true,false,false;kwargs...)
init_SPDE_amimc(;kwargs...) = init_SPDE(AD(2),false,false,false,false;kwargs...)
init_SPDE_amimc_multiple(;kwargs...) = init_SPDE(AD(2),false,true,false,false;kwargs...)
init_SPDE_amiqmc(;kwargs...) = init_SPDE(AD(2),true,false,false,false;kwargs...)
init_SPDE_amiqmc_multiple(;kwargs...) = init_SPDE(AD(2),true,true,false,false;kwargs...)
init_SPDE_mgmlmc(;kwargs...) = init_SPDE(MGML(),false,false,false,true;kwargs...)
init_SPDE_mgmlmc_multiple(;kwargs...) = init_SPDE(MGML(),false,true,false,true;kwargs...)

## main function ##
function init_SPDE(method::IndexSet, is_qmc::Bool, is_multiple_qoi::Bool, is_analyse::Bool, is_multigrid::Bool; corr_len::T=0.5, smoothness::T=1.5, nterms::N=500, max_level::N=5, continuate::Bool=false, nshifts::N=20, verbose::Bool=true, sample_multiplication_factor::T=1.1) where {T<:AbstractFloat,N<:Integer}

    ## Gaussian random fields ##
    nlevels = isa(method,SL) ? 1 : max_level + 1
    coarse_dof = 2

    # covariance function
    cov_fun = CovarianceFunction(2,Matern(corr_len,smoothness))

    # level 0
    m = coarse_dof*2^(max_level-nlevels+1)
    v = 1/2/m:1/2/m:1-1/2/m
    grf = GaussianRandomField(cov_fun,KarhunenLoeve(nterms),v,v)

    # for AD index set, pick all fields in TD manner
    nmethod = isa(method,AD) ? TD(2) : method

    # create fields
    zero_idx = get_index_set(nmethod,0)[1]
    fields = Dict{typeof(zero_idx),typeof(grf)}()
    fields[zero_idx] = grf

    # all other levels
    for idx in get_index_set(nmethod,nlevels-1)
        if !haskey(fields,idx) # avoid duplication of zero_idx
            i = idx[1]
            j = length(idx) > 1 ? idx[2] : i
            m = coarse_dof*2^i
            n = coarse_dof*2^j
            vx = 1/2/m:1/2/m:1-1/2/m
            vy = 1/2/n:1/2/n:1-1/2/n
            grf = GaussianRandomField(cov_fun,KarhunenLoeve(nterms),vx,vy)
            fields[idx] = grf
        end
    end

    # user data
    user_data = SPDE_Data(fields)

    # number generator
    if is_qmc
        number_generator = NormalQMCGenerator(nterms,nshifts)
    else
        number_generator = NormalMCGenerator(nterms)
    end

    # name
    name = "SPDE "
    name = is_analyse ? string(name,"analyse ") : name
    name = isa(method,AD) ? string(name,"A") : name
    name = isa(method,MGML) ? string(name,"MG-") : name
    name = isa(method,ML) || isa(method,MGML) ? string(name,"ML") : MultilevelEstimators.ndims(method) > 1 ? string(name,"MI") : name
    name = is_qmc ? string(name,"Q") : name
    name = string(name,"MC")
    name = is_multiple_qoi ? string(name," (multiple)") : name

    # sample function
    sample_function = "SPDE"
    sample_function = is_multigrid ? string(sample_function,"_mg") : sample_function
    sample_function = is_multiple_qoi ? string(sample_function,"_multiple") : string(sample_function,"_single")
    sample_function = eval(:($(Symbol(sample_function))))

    ## Estimator ##
    create_estimator(
        name = name, # estimator name
        folder = string(joinpath(Pkg.dir("MultilevelEstimators"),"applications","SPDE","data",name)), # for report
        method = method, # method: ML, SL, TD...
        number_generator = number_generator, # number generator
        sample_function = sample_function, # qoi
        user_data = user_data, # GRF's
        verbose = true, # display information
        max_level = nlevels-1, # maximum number of levels
        continuate = continuate, # continuate on larger tolerances
        nb_of_qoi = is_multiple_qoi ? 20^2 : 1, # number of qoi
        cost_model = (index) -> geometric_cost_model(4,1.5,index), # cost model
        sample_multiplication_factor = sample_multiplication_factor # qmc multiplication factor
    )
end
