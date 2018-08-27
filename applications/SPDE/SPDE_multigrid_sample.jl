## SPDE_multigrid_sample.jl : sample functions for lognormal diffusion problem, with multigrid

## solver ##
function SPDE_mg_solve(Z::Matrix{T}) where {T<:Real}
    A = elliptic2d(exp.(Z))
    b = ones(size(A,1)) # rhs
    mg = V_cycle(A,size(Z)) # mg structure
    mg.grids[1].b .= b # copy rhs
    FMG(mg.grids,2,2,1,1,FullWeighting()), getfield.(mg.grids,:sz)
end

function FMG(grids::Vector{G} where {G<:Grid}, ν₀::Int, ν₁::Int, ν₂::Int, grid_ptr::Int, smoother::Smoother)
    if grid_ptr == length(grids)
        grids[grid_ptr].x .= zeros(grids[grid_ptr].x)
        sol = Vector{Vector{eltype[grids[1].x]}}()
    else
        grids[grid_ptr+1].b .= grids[grid_ptr].R*grids[grid_ptr].b
        sol = FMG(grids,ν₀,ν₁,ν₂,grid_ptr+1,smoother,sol)
        psuh!(sol,grids[grid_ptr+1].x) 
        grids[grid_ptr].x .= P(Cubic(),grids[grid_ptr+1].sz...)*grids[grid_ptr+1].x # FMG with cubic interpolation
    end
    for i in 1:ν₀
        μ_cycle!(grids,1,ν₁,ν₂,grid_ptr,smoother)
    end
    return sol
end

## sample functions ##
for mode in ["single" "multiple"]
    ex = :(
           function $(Symbol("lognormal_diffusion_mg_",mode))(index::Index, ξ::Vector{T} where {T<:Real}, data::SPDE_Data)

               # extract grf
               grf = data[index]

               # solve
               Zf = sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
               Qs = $(Symbol("SPDE_mg_sample_",mode))(Zf)

               # compute difference
               for (key,value) in diff(index)
                   Zc = interpolate_field(data[index].pts,data[key].pts,Zf) # interpolation of fine grid GRF
                   Qc = $(Symbol("SPDE_mg_sample_",mode))(Zc)
                   dQ += value*Qc
               end

               return (dQ,Qf)
           end
          )
    eval(ex)
    ex = :(
           function $(Symbol("SPDE_mg_sample_",mode))(Z::Matrix{T}) where {T<:Real}
               (sol,szs) = $(Symbol("SPDE_mg_solve"))(Z::Matrix{T}) where {T<:Real}
               $(Symbol("SPDE_",mode,"_qoi")).(sol,szs)
           end
          )
    eval(ex)
end






