 module Solver_L

using SparseArrays
using LinearAlgebra


function Interface(E_in::Matrix{T},nelx::Int64,nely::Int64,nu_in::T) where {T<:Real}



#println(nelx)
#println(nely)


#println(typeof(nelx))
ndof=2*(nely+1)*(nelx+1)

ndof=Int(ndof)

fixeddofs = 1:1:2*(nely+1);
#println(fixeddofs)
fixeddofsright=2*((nelx+1)*(nely+1))-nely*2-1:1:(nelx+1)*(nely+1)*2;
#println(fixeddofsright)
fixeddofs=[collect(fixeddofs);collect(fixeddofsright)];

Force=1e7;


Fdof = ((nely+1)*nelx)*ones(length(collect(2*(1:nely+1))))+collect(2*(1:nely+1));
Fval = -Force/nely.*[1/2; ones(nely-1); 1/2];
#println(typeof(Fdof))


#println(typeof(ndof))
#println(ndof)
#println(Fdof)

#Fdof=[202;204;206;208;210]
#println(size(Fdof))
sz=size(Fdof)
One=vec(ones(sz[1],1))


F = sparse(collect(Fdof),One,Fval,ndof,1); #matlab
#println("Succes")
#println(F)



xN = ones(nely,nelx);
nu = nu_in*xN;



if(size(E_in,2)==1)
E=E_in*xN;
else
E=Array(transpose(E_in));
end
#println(typeof(fixeddofs))

u = planestress4(E,nu,fixeddofs,F)

return u

end

function planestress4(E::Matrix{T},nu::Matrix{T},fixeddofs::Array{Int64},F::SparseMatrixCSC{T}) where {T<:Real}

#     println("started")
    # Number of elements
    nelx = Int(size(E,2));
    nely = Int(size(E,1));


    # Degrees of freedom
    ndof = 2*(nely+1)*(nelx+1);
    alldofs = 1:ndof;
    freedofs = setdiff(alldofs,fixeddofs);
#    println("1")

    # Matrices with degrees of freedom and nodes for each element
    nodenrs = reshape(1:(nelx+1)*(nely+1),nely+1,nelx+1);
    intermNodenrs=nodenrs[1:end-1,(1:end-1)].+1
#    println("2")
#    println(typeof(intermNodenrs))
    mult=nelx*nely

    edofVec = reshape((2*intermNodenrs),mult,1);
#    println("3")

    edofMat = repeat(edofVec,1,8).+repeat([0 1 2*nely.+[2 3 0 1] -2 -1],nelx*nely,1);
#    println("4")

    enodeVec = reshape(nodenrs[2:end,1:end-1],nelx*nely,1);
#    println("5")

    enodeMat = repeat(enodeVec,1,4).+repeat([0 nely.+[1 0] -1],nelx*nely,1);
#    println("6")


    # Number of (neighbour) elements for each node.
    neighb = 4*ones(1+nely,1+nelx);
    neighb[1,:] = neighb[1,:]/2;
    neighb[end,:] = neighb[end,:]/2;
    neighb[:,1] = neighb[:,1]/2;
    neighb[:,end] = neighb[:,end]/2;
    neighb = neighb[:];
    nele = neighb[enodeMat];

#    println("10")
    # Stiffness matrix
#println(typeof(kron(edofMat,ones(8))'))

iK = reshape(kron(edofMat,ones(8))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(8)')',64*nelx*nely,1);

A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
A = [A11 A12; transpose(A12) A11];
B = [B11 B12; transpose(B12) B11];
#println("11")


sKe = A[:]*(Array(transpose(E[:]))./(1 .- Array(transpose(nu[:])).^2)/24).+B[:]*(Array(transpose(nu[:])).*Array(transpose(E[:]))./(1 .- Array(transpose(nu[:])).^2)/24);
#sKe=1
sK = reshape(sKe,64*nelx*nely,1);

#println(sK[1161])
#println(size(jK))
#println(size(iK))


#iK=Int(iK)
#jk=Int(jK)

#println(jK)

iK=iK.-1
jK=jK.-1


K = sparse(vec(iK),vec(jK),vec(sK));
#if(condest(K)<1.5e-17)
#    d=0;
#end

#println("12")


# Displacements
u = zeros(ndof,1);
K_int=K[freedofs,freedofs]
F_int=F[freedofs,1]
#println(size(K_int))
#println(size(F_int))

K=0


f(K_int,F_int)=try
 K_int\vec(Array(F_int)) ;
catch
    println("Error Thrown")
    println("Compute with different method")
end



u_int=f(K_int,F_int)




u[freedofs,1]=u_int

F_int=0
K_int=0
u_int=0
GC.gc()
#println("13")
#println(u)
#F_r=full(F_int)
#println(F_r)
#K_r=full(K)
#println(K_r[18,9])

#println(u[210])
#println(K)


return u

end
end
