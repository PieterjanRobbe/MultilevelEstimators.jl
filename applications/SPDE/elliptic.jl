## elliptic.jl : finite difference discretization of general elliptic PDE in 1D and 2D

function estencil(k)
    d1 = k[1:end-1]+k[2:end]
    d2 = -k[2:end-1]
    return spdiagm((d2,d1,d2),(-1,0,1))
end

"""
elliptic1d(k)

Generate system matrix for general 1d elliptic PDE where k is the value of the diffusion coefficient at the mid points
"""
function elliptic1d(k)
    return length(k)^2*estencil(k)
end

"""
elliptic2d(kx,ky)

Generate system matrix for general 2d elliptic PDE on an m-by-n grid where k is the value of the diffusion coefficient at the mid points of the grid.

kx is an m-by-(n-1) array that contains the k_(i+1/2,j) values
ky is an (m-1)-by-n array that contains the k_(i,j+1/2) values
"""
function elliptic2d(kx,ky)
    m = size(kx,1)
    n = size(ky,2)
    B1 = blkdiag([m^2*estencil(kx[:,j]) for j=1:n-1]...)
    B2 = blkdiag([n^2*spdiagm(ky[:,i]+ky[:,i+1]) for i=1:n-1]...)
    c = -vcat([n^2*ky[:,i] for i=2:n-1]...)	
    C1 = spdiagm((c,c),(m-1,-(m-1)))	
    return isempty(C1) ? B1+B2 : B1+B2+C1
end

"""
elliptic2d(k)

Generate system matrix for general 2d elliptic PDE on an m-by-n grid where k is the value of the diffusion coefficient at the mid points of the grid.

k is an m-by-n array that contains the k values (of which only half will be used)
"""
elliptic2d(k) = elliptic2d(k[1:2:end,2:2:end-1],k[2:2:end-1,1:2:end])
