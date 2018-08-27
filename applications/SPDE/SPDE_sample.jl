## SPDE_sample.jl : sample functions for lognormal diffusion problem

## solver ##
function SPDE_solve(Z::Matrix{T}) where {T<:Real}
    A = elliptic2d(exp.(Z))
    b = ones(size(A,1))
    A\b
end

## qoi ##
function SPDE_single_qoi(x::Vector{T},sz::Tuple{Int,Int}) where {T<:Real}
    m,n = round.(Int,(sz.+1)./2)
    i = sub2ind((m-1,n-1),ceil.(Int,(m,n)./2)...)
    return x[i]
end

function SPDE_multiple_qoi(x::Vector{T},sz::Tuple{Int,Int}) where {T<:Real}
    m,n = round.(Int,(sz.+1)./2).-1
    x_reshaped = reshape(x,(m,n))
    x_padded = hcat(zeros(m,1),x_reshaped,zeros(m,1)) # pad solution with dirichlet conditions
    x_padded = vcat(zeros(1,n+2),x_padded,zeros(1,n+2))
    itp = interpolate(linspace.(0,1,(m+2,n+2)), x_padded, Gridded(Linear()))
    pts = linspace(0,1,20)
    return itp[pts,pts][:]
end

## interpolate field ##
function interpolate_field(pts_fine,pts_coarse,Z::Matrix{T}) where {T<:Real}
    itp = interpolate(pts_fine, Z, Gridded(Linear()))
    itp[pts_coarse[1],pts_coarse[2]]
end

## sample functions ##
for mode in ["single" "multiple"]
    ex = :(
           function $(Symbol("lognormal_diffusion_",mode))(index::Index, ξ::Vector{T} where {T<:Real}, data::SPDE_Data)

               # extract grf
               grf = data[index]

               # solve
               Zf = sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
               Qf = $(Symbol("SPDE_sample_",mode))(Zf)

               # compute difference
               dQ = Qf
               for (key,value) in diff(index)
                   Zc = interpolate_field(data[index].pts,data[key].pts,Zf) # interpolation of fine grid GRF
                   Qc = $(Symbol("SPDE_sample_",mode))(Zc)
                   dQ += value*Qc
               end

               return (dQ,Qf)
           end
          )
    eval(ex)
    ex = :( $(Symbol("SPDE_sample_",mode))(Z::Matrix{T}) where {T<:Real} = $(Symbol("SPDE_",mode,"_qoi"))(SPDE_solve(Z),size(Z)) )
    eval(ex)
end
