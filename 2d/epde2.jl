######################################################################################################
## solve elliptic PDE in two dimensions using a FV approach
######################################################################################################
# The PDE must be of the form
#        -d/dx(k(x)dp/dx) = f(x) x \in [0,1]^2
# where p is unknown, f is a RHS and k(x) is the diffusion coefficient.
# Boundary conditions are of Dirichlet type at x = 0 and x = 1, and of
# Neumann type at y = 0 and y = 1. Here, the boundary
# conditions are 1 at x = 0 and 0 everywhere else. There is no source term (f = 0).
# INPUT:
#   k = array containing the value of k at the cell centers x_(i,j), (i,j) = 1..m
# OUTPUT:
#   matrix containing the value of p at the cell centers x_(i,j), (i,j) = 1..m
function epde2{T<:Number}(k::Array{T,2})
  # number of meshes
  mx = size(k,1)
  my = size(k,2)

  # harmonic mean values in x-direction
  k_hx = zeros(k)
  dx = mx/my
  for i = 1:mx-1
    k_hx[i,:,:] = harmmean(k[i:i+1,:],1)*dx
  end

  # harmonic mean values in x-direction
  k_hy = zeros(k)
  dy = my/mx
  for i = 1:my-1
    k_hy[:,i,:] = harmmean(k[:,i:i+1],2)*dy
  end

  # main diagonal
  sigma = zeros(mx,my);

  # bottom layer ( k = 1 )
  sigma[1,1]           = k_hx[1,1]           + 2*k[1,1]*dx         + k_hy[1,1]           + 0
  sigma[2:mx-1,1]      = k_hx[2:mx-1,1]      + k_hx[1:mx-2,1]      + k_hy[2:mx-1,1]      + 0
  sigma[mx,1]          = k_hx[mx-1,1]        + 2*k[mx,1]*dx        + k_hy[mx,1]          + 0
  sigma[1,2:my-1]      = k_hx[1,2:my-1]      + 2*k[1,2:my-1]*dx    + k_hy[1,2:my-1]      + k_hy[1,1:my-2]
  sigma[2:mx-1,2:my-1] = k_hx[2:mx-1,2:my-1] + k_hx[1:mx-2,2:my-1] + k_hy[2:mx-1,2:my-1] + k_hy[2:mx-1,1:my-2]
  sigma[mx,2:my-1]     = k_hx[mx-1,2:my-1]   + 2*k[mx,2:my-1]*dx   + k_hy[mx,2:my-1]     + k_hy[mx,1:my-2]
  sigma[1,my]          = k_hx[1,my]          + 2*k[1,my]*dx        + k_hy[1,my-1]        + 0
  sigma[2:mx-1,my]     = k_hx[2:mx-1,my]     + k_hx[1:mx-2,my]     + k_hy[2:mx-1,my-1]   + 0
  sigma[mx,my]         = k_hx[mx-1,my]       + 2*k[mx,my]*dx       + k_hy[mx,my-1]       + 0

  # transform into vector
  sigma = sigma[:]

  # off-diagonal elements
  k_hx[mx,:] = 0.
  k_hy[:,my] = 0.

  # compose sparse matrix
  F = spdiagm((sigma,-k_hx[1:mx*my-1],-k_hx[1:mx*my-1],-k_hy[1:mx*my-mx],-k_hy[1:mx*my-mx]),(0,1,-1,mx,-mx))

  # compose right-hand side
  b = zeros(mx,my)
  b[1,:] = 2*k[1,:]*dx
  b = b[:]

  # solve
  return reshape(F\b[:],mx,my)
end

######################################################################################################
## solve elliptic PDE in two dimensions using a FV approach
######################################################################################################
# The PDE must be of the form
#        -d/dx(k(x)dp/dx) = f(x) x \in [0,1]^3
# where p is unknown, f is a RHS and k(x) is the diffusion coefficient.
# Boundary conditions are of Dirichlet type at x = 0, x = 1, y = 0 and y = 1.
# Here, the boundary conditions are 1 at x = 0 and 0 everywhere else. There is no source term (f = 0).
# INPUT:
#   k = array containing the value of k at the cell centers x_(i,j), (i,j) = 1..m
# OUTPUT:
#   matrix containing the value of p at the cell centers x_(i,j), (i,j) = 1..m
function epde2d{T<:Number}(k::Array{T,2})
  # number of meshes
  mx = size(k,1)
  my = size(k,2)

  # harmonic mean values in x-direction
  k_hx = zeros(k)
  dx = mx/my
  for i = 1:mx-1
    k_hx[i,:,:] = harmmean(k[i:i+1,:],1)*dx
  end

  # harmonic mean values in x-direction
  k_hy = zeros(k)
  dy = my/mx
  for i = 1:my-1
    k_hy[:,i,:] = harmmean(k[:,i:i+1],2)*dy
  end

  # main diagonal
  sigma = zeros(mx,my);

  # bottom layer ( k = 1 )
  sigma[1,1]           = k_hx[1,1]           + 2*k[1,1]*dx         + k_hy[1,1]           + 2*k[1,1]*dy
  sigma[2:mx-1,1]      = k_hx[2:mx-1,1]      + k_hx[1:mx-2,1]      + k_hy[2:mx-1,1]      + 2*k[2:mx-1,1]*dy 
  sigma[mx,1]          = k_hx[mx-1,1]        + 2*k[mx,1]*dx        + k_hy[mx,1]          + 2*k[mx,1]*dy

  sigma[1,2:my-1]      = k_hx[1,2:my-1]      + 2*k[1,2:my-1]*dx    + k_hy[1,2:my-1]      + k_hy[1,1:my-2]
  sigma[2:mx-1,2:my-1] = k_hx[2:mx-1,2:my-1] + k_hx[1:mx-2,2:my-1] + k_hy[2:mx-1,2:my-1] + k_hy[2:mx-1,1:my-2]
  sigma[mx,2:my-1]     = k_hx[mx-1,2:my-1]   + 2*k[mx,2:my-1]*dx   + k_hy[mx,2:my-1]     + k_hy[mx,1:my-2]

  sigma[1,my]          = k_hx[1,my]          + 2*k[1,my]*dx        + k_hy[1,my-1]        + 2*k[1,my]*dy
  sigma[2:mx-1,my]     = k_hx[2:mx-1,my]     + k_hx[1:mx-2,my]     + k_hy[2:mx-1,my-1]   + 2*k[2:mx-1,my]*dy
  sigma[mx,my]         = k_hx[mx-1,my]       + 2*k[mx,my]*dx       + k_hy[mx,my-1]       + 2*k[mx,my]*dy

  # transform into vector
  sigma = sigma[:]

  # off-diagonal elements
  k_hx[mx,:] = 0.
  k_hy[:,my] = 0.

  # compose sparse matrix
  F = spdiagm((sigma,-k_hx[1:mx*my-1],-k_hx[1:mx*my-1],-k_hy[1:mx*my-mx],-k_hy[1:mx*my-mx]),(0,1,-1,mx,-mx))

  # compose right-hand side
  b = ones(mx,my)/mx/my
  b = b[:]

  # solve
  return reshape(F\b[:],mx,my)
end

# compute the harmonic mean of A along dimension d
function harmmean{T<:Number,S<:Integer}(A::AbstractArray{T}, d::S)
    return squeeze(1 ./ mean( 1./A, d), d)
end

