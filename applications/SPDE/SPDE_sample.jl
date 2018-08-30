## SPDE_sample.jl : sample functions for lognormal diffusion problem

## solver ##
function SPDE_solve(Z::Matrix{T}) where {T<:Real}
	A = elliptic2d(exp.(Z))
    b = ones(size(A,1))
    A\b
end

## qoi ##
SPDE_single_qoi(x::Vector{T},sz::Tuple{Int,Int}) where {T<:Real} = x[round(Int,(prod(sz.-1)+1)/2)]

function SPDE_multiple_qoi(x::Vector{T},sz::Tuple{Int,Int}) where {T<:Real}
    x_reshaped = reshape(x,sz.-1)
	x_padded = PaddedView(0, x_reshaped, sz.+1, (2,2))
    itp = interpolate(linspace.(0,1,sz.+1), x_padded, Gridded(Linear()))
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
           function $(Symbol("SPDE_sample_",mode))(index::Index, ξ::Vector{T} where {T<:Real}, data::SPDE_Data)

               # extract grf
               grf = data[index]

               # solve
               Zf = sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
               Qf = $(Symbol("SPDE_single_sample_",mode))(Zf)

               # compute difference
               dQ = Qf
               for (key,value) in diff(index)
                   Zc = interpolate_field(data[index].pts,data[key].pts,Zf) # interpolation of fine grid GRF
                   Qc = $(Symbol("SPDE_single_sample_",mode))(Zc)
                   dQ += value*Qc
               end

			   # safety
			   while !is_valid_sample(Qf)
			   	   	 ξ₂ = randn(size(ξ))
               		 Zf = sample(grf,xi=ξ₂[1:randdim(grf)]) # compute GRF
               		 Qf = $(Symbol("SPDE_single_sample_",mode))(Zf)
					 dQ = Qf
                	 for (key,value) in diff(index)
                   		 Zc = interpolate_field(data[index].pts,data[key].pts,Zf) # interpolation of fine grid GRF
                   		 Qc = $(Symbol("SPDE_single_sample_",mode))(Zc)
                   		 dQ += value*Qc
               		 end
			   end

               return (dQ,Qf)
           end
          )
    eval(ex)
    ex = :( $(Symbol("SPDE_single_sample_",mode))(Z::Matrix{T}) where {T<:Real} = $(Symbol("SPDE_",mode,"_qoi"))(SPDE_solve(Z),size(Z).-1) )
    eval(ex)
end

is_valid_sample(Qf) = !(any(j->j>1||j<0,Qf))
