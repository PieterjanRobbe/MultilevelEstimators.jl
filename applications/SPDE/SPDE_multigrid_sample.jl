## SPDE_multigrid_sample.jl : sample functions for lognormal diffusion problem, with multigrid

## solver ##
function SPDE_mg_solve(Z::Matrix{T},index::Index) where {T<:Real}
    A = elliptic2d(exp.(Z))

    b = ones(size(A,1)) # rhs
    if length(index) > 1
        mg = NotSoSimpleMultigrid.V_cycle(A,size(Z).-1) # mg structure
    else
        mg = SimpleMultigrid.V_cycle(A,size(Z).-1) # mg structure
    end
    mg.grids[1].b .= b # copy rhs

	FMG(mg.grids,2,2,1,1,GaussSeidel()), getfield.(mg.grids,:sz)[StepRange.(size(mg.grids),-1,size(mg.grids))...]
end

## sample functions ##
for mode in ["single" "multiple"]
    ex = :(
           function $(Symbol("SPDE_sample_mg_",mode))(index::Index, ξ::Vector{T} where {T<:Real}, data::SPDE_Data)

               # extract grf
               grf = data[index]

               # solve
               Zf = sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
               Qf = $(Symbol("SPDE_single_sample_mg_",mode))(Zf,index)

               # safety
               while !(is_valid_sample_mg(Qf))
                   Zf = sample(grf,xi=randn(size(ξ))) # recompute GRF
                   Qf = $(Symbol("SPDE_single_sample_mg_",mode))(Zf,index)
               end

               # compute difference
               dQ = deepcopy(Qf)
			   for ℓ in Iterators.product(colon.(0,size(dQ).-1)...)
                   index_ = Index(ℓ...)
                   for (key,value) in diff(index_)
                       dQ[(index_.+1)...] += value*Qf[(key.+1)...]
                   end
               end

               return (dQ,Qf)
           end
          )
    eval(ex)
    ex = :(
           function $(Symbol("SPDE_single_sample_mg_",mode))(Z::Matrix{T},index::Index) where {T<:Real}
               (sol,szs) = $(Symbol("SPDE_mg_solve"))(Z,index)
               $(Symbol("SPDE_",mode,"_qoi")).(sol,szs)
           end
          )
    eval(ex)
end

function is_valid_sample_mg(Qf)
    check = true
    for i in 1:length(Qf)
        if any(j->j>1||j<0,Qf[i])
            check = false
        end
    end 
    return check
end
