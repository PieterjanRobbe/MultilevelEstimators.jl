# solver.jl : solver for 2d lognormal diffusion problem

function lognormal_diffusion(index::Index, ξ::Vector{T} where {T<:Real}, sampler::Sampler)
    # extract grf
    grf = sampler.user_data[index]

    # solve
    A = elliptic2d(exp.(sample(grf,xi=ξ)))
    b = ones(size(A,1))
    x = A\b

    # compute qoi
    m,n = round.(Int,(length.(grf.pts).+1)./2)
    i = sub2ind((m-1,n-1),ceil.(Int,(m,n)./2)) 

    return x[i]
end

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

# init functions
init_lognormal_diffusion_mc() = init_lognormal_diffusion(0,false)

function init_lognormal_diffusion(method::MLE.IndexSet, is_qmc::Bool)

    ## Gaussian random fields ##
    corr_len = 0.2
    smoothness = 2.0
    nterms = 1000
    max_level = isa(method,SL) ? 1 : 6

    # covariance function
    cov_fun = CovarianceFunction(2,Matern(corr_len,smoothness))

    # level 0
    m = 2*2^max_level
    v = 1/2/m:1/2/m:1-1/2/m
    grf = GaussianRandomField(cov_fun,KarhunenLoeve(nterms),v,v)

    # create fields
    zero_idx = get_index_set(method,0)[1]
    fields = Dict{typeof(zero_idx),typeof(grf)}()
    fields[zero_idx] = grf

    # all other levels
    for idx in get_index_set(method,max_level)
        if !haskey(fields,idx) # avoid duplication of zero_idx
            i = idx[1]
            j = length(idx) > 1 ? idx[2] : i
            m = 2*2^i
            n = 2*2^j
            vx = 1/2/m:1/m:1-1/2/m
            vy = 1/2/n:1/n:1-1/2/n
            grf = GaussianRandomField(cov_fun,KarhunenLoeve(nterms),vx,vy)
            fields[idx] = grf
        end
    end

    # user data
    user_data = SPDE_Data(fields)

    ## Estimator ##
    ##
    ##
    create_estimator(method=method,number_generator=UniformMCGenerator(100),sample_function=(i)->randn(i),tol=1e-2)
    ##
    ##

end

struct SPDE_Data where V
    fields::V
end
