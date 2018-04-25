module SPDE

using Interpolations, Reexport, PyPlot
@reexport using MultilevelEstimators, GaussianRandomFields

import Base.getindex

export init_lognormal_diffusion_analyse_ml
export init_lognormal_diffusion_mc
export init_lognormal_diffusion_mc_multiple
export init_lognormal_diffusion_mlmc
export init_lognormal_diffusion_mlmc_multiple
export init_lognormal_diffusion_qmc
export init_lognormal_diffusion_qmc_multiple
export init_lognormal_diffusion_mlqmc
export init_lognormal_diffusion_mlqmc_multiple

# user data type to hold GRF's
struct SPDE_Data{V}
    fields::V
end

getindex(s::SPDE_Data,index::Index) = s.fields[index]

## init functions ##

init_lognormal_diffusion_analyse_ml() = init_lognormal_diffusion(ML(),false,false,true)
init_lognormal_diffusion_mc() = init_lognormal_diffusion(SL(),false,false,false)
init_lognormal_diffusion_mc_multiple() = init_lognormal_diffusion(SL(),false,true,false)
init_lognormal_diffusion_mlmc() = init_lognormal_diffusion(ML(),false,false,false)
init_lognormal_diffusion_mlmc_multiple() = init_lognormal_diffusion(ML(),false,true,false)
init_lognormal_diffusion_qmc() = init_lognormal_diffusion(SL(),true,false,false)
init_lognormal_diffusion_qmc_multiple() = init_lognormal_diffusion(SL(),true,true,false)
init_lognormal_diffusion_mlqmc() = init_lognormal_diffusion(ML(),true,false,false)
init_lognormal_diffusion_mlqmc_multiple() = init_lognormal_diffusion(ML(),true,true,false)

function init_lognormal_diffusion(method::IndexSet, is_qmc::Bool, is_multiple_qoi::Bool, is_analyse::Bool)

    ## Gaussian random fields ##
    corr_len = 0.5
    smoothness = 1.5
    nterms = 500
    max_level = 4
    nlevels = isa(method,SL) ? 1 : max_level + 1
    coarse_dof = 2

    # covariance function
    cov_fun = CovarianceFunction(2,Matern(corr_len,smoothness))

    # level 0
    m = coarse_dof*2^(max_level-nlevels+1)
    v = 1/2/m:1/2/m:1-1/2/m
    grf = GaussianRandomField(cov_fun,KarhunenLoeve(nterms),v,v)

    # create fields
    zero_idx = get_index_set(method,0)[1]
    fields = Dict{typeof(zero_idx),typeof(grf)}()
    fields[zero_idx] = grf

    # all other levels
    for idx in get_index_set(method,nlevels-1)
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
        nshifts = 20
        number_generator = NormalQMCGenerator(nterms,nshifts)
    else
        number_generator = NormalMCGenerator(nterms)
    end

    # name
    name = "SPDE "
    name = is_analyse ? string(name,"analyse ") : name
    name = isa(method,ML) ? string(name,"ML") : name
    name = is_qmc ? string(name,"Q") : name
    name = string(name,"MC")
    name = is_multiple_qoi ? string(name," (multiple)") : name

    ## Estimator ##
    create_estimator(
        name = name,
        folder = string(Pkg.dir("MultilevelEstimators"),"applications/SPDE/data/",name),
        method = method,
        number_generator = number_generator,
        sample_function = is_multiple_qoi ? lognormal_diffusion_multiple : lognormal_diffusion_single,
        user_data = user_data,
        verbose = true,
        max_level = nlevels-1,
        continuate = true,
        nb_of_qoi = is_multiple_qoi ? 20^2 : 1,
        #do_regression=false,
        cost_model = (index) -> geometric_cost_model(4,1.5,index),
        #do_splitting = false,
        #sample_multiplication_factor = 1.1
    )
end


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
    x_padded = vcat(zeros(1,m+2),x_padded,zeros(1,m+2))
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

    # points where to compute the solution
    max_idx = sort(collect(keys(data.fields)))[end] # find largest grid
    pts_ = data.fields[max_idx].pts
    m = round.(Int,(length.(pts_).+1)./2).+1

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
