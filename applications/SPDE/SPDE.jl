#=
SPDE.jl : Module for simulating a PDE with random coefficents using various
          methods as implemented in MultilevelEstimators.jl.
This file defines a module that contains methods to initialize the estimator.
Currently implemented:
    - MC        : Monte Carlo
    - QMC       : Quasi-Monte Carlo
    - MLMC      : Multilevel Monte Carlo
    - MLQMC     : Multilevel Quasi-Monte Carlo
    - MIMC      : Multi-Index Monte Carlo
    - MIQMC     : Multi-Index Quasi-Monte Carlo
    - AMIMC     : Adaptive Multi-Index Monte Carlo
    - AMIQMC    : Adaptive Multi-Index Quasi-Monte Carlo
NOTE: Technically, we are simulating a PDE with random coefficients, but this is too
long to write down, so with a slight abuse of notation, this module is called "SPDE".
=#
module SPDE

## dependencies ##
using Interpolations, Reexport, PyPlot
@reexport using MultilevelEstimators, GaussianRandomFields

## import statements ##
import Base.getindex

## export statements ##
export init_lognormal_diffusion_analyse_ml
export init_lognormal_diffusion_analyse_mi
export init_lognormal_diffusion_mc
export init_lognormal_diffusion_mc_multiple
export init_lognormal_diffusion_mlmc
export init_lognormal_diffusion_mlmc_multiple
export init_lognormal_diffusion_qmc
export init_lognormal_diffusion_qmc_multiple
export init_lognormal_diffusion_mlqmc
export init_lognormal_diffusion_mlqmc_multiple
export init_lognormal_diffusion_mimc
export init_lognormal_diffusion_mimc_multiple
export init_lognormal_diffusion_miqmc
export init_lognormal_diffusion_miqmc_multiple
export init_lognormal_diffusion_amimc
export init_lognormal_diffusion_amimc_multiple
export init_lognormal_diffusion_amiqmc
export init_lognormal_diffusion_amiqmc_multiple

## init functions ##
init_lognormal_diffusion_analyse_ml(;kwargs...) = init_lognormal_diffusion(ML(),false,false,true;kwargs...)
init_lognormal_diffusion_analyse_mi(;kwargs...) = init_lognormal_diffusion(TD(2),false,false,true;kwargs...)
init_lognormal_diffusion_mc(;kwargs...) = init_lognormal_diffusion(SL(),false,false,false;kwargs...)
init_lognormal_diffusion_mc_multiple(;kwargs...) = init_lognormal_diffusion(SL(),false,true,false;kwargs...)
init_lognormal_diffusion_mlmc(;kwargs...) = init_lognormal_diffusion(ML(),false,false,false;kwargs...)
init_lognormal_diffusion_mlmc_multiple(;kwargs...) = init_lognormal_diffusion(ML(),false,true,false;kwargs...)
init_lognormal_diffusion_qmc(;kwargs...) = init_lognormal_diffusion(SL(),true,false,false;kwargs...)
init_lognormal_diffusion_qmc_multiple(;kwargs...) = init_lognormal_diffusion(SL(),true,true,false;kwargs...)
init_lognormal_diffusion_mlqmc(;kwargs...) = init_lognormal_diffusion(ML(),true,false,false;kwargs...)
init_lognormal_diffusion_mlqmc_multiple(;kwargs...) = init_lognormal_diffusion(ML(),true,true,false;kwargs...)
init_lognormal_diffusion_mimc(;kwargs...) = init_lognormal_diffusion(TD(2),false,false,false;kwargs...)
init_lognormal_diffusion_mimc_multiple(;kwargs...) = init_lognormal_diffusion(TD(2),false,true,false;kwargs...)
init_lognormal_diffusion_miqmc(;kwargs...) = init_lognormal_diffusion(TD(2),true,false,false;kwargs...)
init_lognormal_diffusion_miqmc_multiple(;kwargs...) = init_lognormal_diffusion(TD(2),true,true,false;kwargs...)
init_lognormal_diffusion_amimc(;kwargs...) = init_lognormal_diffusion(AD(2),false,false,false;kwargs...)
init_lognormal_diffusion_amimc_multiple(;kwargs...) = init_lognormal_diffusion(AD(2),false,true,false;kwargs...)
init_lognormal_diffusion_amiqmc(;kwargs...) = init_lognormal_diffusion(AD(2),true,false,false;kwargs...)
init_lognormal_diffusion_amiqmc_multiple(;kwargs...) = init_lognormal_diffusion(AD(2),true,true,false;kwargs...)

function init_lognormal_diffusion(method::IndexSet, is_qmc::Bool, is_multiple_qoi::Bool, is_analyse::Bool; corr_len::T=0.5, smoothness::T=1.5, nterms::N=500, max_level::N=5, continuate::Bool=false, nshifts::N=20, verbose::Bool=true, sample_multiplication_factor::T=1.1) where {T<:AbstractFloat,N<:Integer}

    ## Gaussian random fields ##
    nlevels = isa(method,SL) ? 1 : max_level + 1
    @show nlevels
    @show method
    @show SL
    @show isa(method,SL)
    coarse_dof = 2

    # covariance function
    cov_fun = CovarianceFunction(2,Matern(corr_len,smoothness))

    # level 0
    m = coarse_dof*2^(max_level-nlevels+1)
    @show m
    v = 1/2/m:1/2/m:1-1/2/m
    @show v
    grf = GaussianRandomField(cov_fun,KarhunenLoeve(nterms),v,v)
    plot(grf)

    # for AD index set, pich all fields in TD manner
    nmethod = isa(method,AD) ? TD(2) : method
    @show nmethod
    @show TD(2)


    # create fields
    zero_idx = get_index_set(nmethod,0)[1]
    @show zero_idx
    fields = Dict{typeof(zero_idx),typeof(grf)}()
    fields[zero_idx] = grf

    # all other levels
    @show  get_index_set(nmethod,nlevels-1)
    for idx in get_index_set(nmethod,nlevels-1)
        @show idx
        if !haskey(fields,idx) # avoid duplication of zero_idx
            i = idx[1]
            j = length(idx) > 1 ? idx[2] : i
            m = coarse_dof*2^i
            n = coarse_dof*2^j
            @show m
            @show n
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
    name = isa(method,ML) ? string(name,"ML") : MultilevelEstimators.ndims(method) > 1 ? string(name,"MI") : name
    name = is_qmc ? string(name,"Q") : name
    name = string(name,"MC")
    name = is_multiple_qoi ? string(name," (multiple)") : name

    ## Estimator ##
    create_estimator(
        name = name, # estimator name
        folder = string(joinpath(Pkg.dir("MultilevelEstimators"),"applications","SPDE","data",name)), # for report
        method = method, # method: ML, SL, TD...
        number_generator = number_generator, # number generator
        sample_function = is_multiple_qoi ? lognormal_diffusion_multiple : lognormal_diffusion_single, # qoi
        user_data = user_data, # GRF's
        verbose = true, # display information
        max_level = nlevels-1, # maximum number of levels
        continuate = continuate, # continuate on larger tolerances
        nb_of_qoi = is_multiple_qoi ? 20^2 : 1, # number of qoi
        cost_model = (index) -> geometric_cost_model(4,1.5,index), # cost model
        sample_multiplication_factor = sample_multiplication_factor # qmc multiplication factor
    )
end

## user data ##
struct SPDE_Data{V}
    fields::V
end

getindex(s::SPDE_Data,index::Index) = s.fields[index]


## sample functions ##
function SPDE_sample(Z::Matrix{T}) where {T<:Real}

    # solve system
    A = elliptic2d(exp.(Z))
    b = ones(size(A,1))
    x = A\b

    # compute qoi
    m,n = round.(Int,(size(Z).+1)./2)
    i = sub2ind((m-1,n-1),ceil.(Int,(m,n)./2)...)
    return x[i]
end


function SPDE_sample_multiple(Z::Matrix{T}) where {T<:Real}

    # solve system
    A = elliptic2d(exp.(Z))
    b = ones(size(A,1))
    x = A\b

    # compute qoi's
    m,n = round.(Int,(size(Z).+1)./2).-1
    x_reshaped = reshape(x,(m,n))
    x_padded = hcat(zeros(m,1),x_reshaped,zeros(m,1)) # pad solution with dirichlet conditions
    x_padded = vcat(zeros(1,n+2),x_padded,zeros(1,n+2))
    itp = interpolate(linspace.(0,1,(m+2,n+2)), x_padded, Gridded(Linear()))
    pts = linspace(0,1,20)
    return itp[pts,pts][:]
end

function interpolate_field(pts_fine,pts_coarse,Z::Matrix{T}) where {T<:Real}
    itp = interpolate(pts_fine, Z, Gridded(Linear()))
    itp[pts_coarse[1],pts_coarse[2]]
end

function lognormal_diffusion_single(index::Index, ξ::Vector{T} where {T<:Real}, data::SPDE_Data)

    # extract grf
    grf = data[index]

    # solve
    Zf = sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
    Qf = SPDE_sample(Zf)

    # compute difference
    dQ = Qf
    for (key,value) in diff(index)
        Zc = interpolate_field(data[index].pts,data[key].pts,Zf) # interpolation of fine grid GRF
        Qc = SPDE_sample(Zc)
        dQ += value*Qc
    end

    return (dQ,Qf)
end

function lognormal_diffusion_multiple(index::Index, ξ::Vector{T} where {T<:Real}, data::SPDE_Data)

    # extract grf
    grf = data[index]

    # solve
    Zf = sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
    Qf = SPDE_sample_multiple(Zf)

    # compute difference
    dQ = Qf
    for (key,value) in diff(index)
        Zc = interpolate_field(data[index].pts,data[key].pts,Zf) # interpolation of fine grid GRF
        Qc = SPDE_sample_multiple(Zc)
        dQ += value*Qc
    end

    return (dQ,Qf)
end

## solver ##

# laplacian
function estencil(k)
    d1 = k[1:end-1]+k[2:end]
    d2 = -k[2:end-1]
    return spdiagm((d2,d1,d2),(-1,0,1))
end

# kx is an m-by-(n-1) array that contains the k_(i+1/2,j) values
# ky is an (m-1)-by-n array that contains the k_(i,j+1/2) values
function elliptic2d(kx,ky)
    m = size(kx,1)
    n = size(ky,2)
    B1 = blkdiag([m^2*estencil(kx[:,j]) for j=1:n-1]...)
    B2 = blkdiag([n^2*spdiagm(ky[:,i]+ky[:,i+1]) for i=1:n-1]...)
    c = -vcat([n^2*ky[:,i] for i=2:n-1]...)
    C1 = spdiagm((c,c),(m-1,-(m-1)))
    return isempty(C1) ? B1+B2 : B1+B2+C1
end

# k is an m-by-n array that contains the k values (of which only three quarters will be used)
elliptic2d(k) = elliptic2d(k[1:2:end,2:2:end-1],k[2:2:end-1,1:2:end])

end # module SPDE
