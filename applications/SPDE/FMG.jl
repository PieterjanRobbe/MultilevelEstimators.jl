# FMG.jl : custom implementation of FMG methods from (NotSo)SimpleMultigrid.jl
# that return the coarse solutions

# MG
function FMG(grids::Vector{G} where {G<:SimpleMultigrid.Grid}, ν₀::Int, ν₁::Int, ν₂::Int, grid_ptr::Int, smoother::SimpleMultigrid.Smoother)
    if grid_ptr == length(grids)
        grids[grid_ptr].x .= zeros(grids[grid_ptr].x)
        sol = Vector{Vector{eltype(grids[1].x)}}(length(grids))
    else
        grids[grid_ptr+1].b .= grids[grid_ptr].R*grids[grid_ptr].b
        sol = FMG(grids,ν₀,ν₁,ν₂,grid_ptr+1,smoother)
        grids[grid_ptr].x .= SimpleMultigrid.P(SimpleMultigrid.Cubic(),grids[grid_ptr+1].sz...)*grids[grid_ptr+1].x # FMG with cubic interpolation
    end
    # V-cycling
    for i in 1:ν₀
        SimpleMultigrid.μ_cycle!(grids,1,ν₁,ν₂,grid_ptr,smoother)
    end
    # safety
    if SimpleMultigrid.norm_of_residu(grids[grid_ptr]) >= 1/prod(grids[grid_ptr].sz)
        grids[grid_ptr].x .= grids[grid_ptr].A\grids[grid_ptr].b # exact solve
    end
    sol[length(grids)-grid_ptr+1] = copy(grids[grid_ptr].x) 
    return sol
end

# MSG
function FMG(grids::Array{G} where {G<:SimpleMultigrid.Grid}, ν₀::Int, ν₁::Int, ν₂::Int, grid_ptr::Int, smoother::SimpleMultigrid.Smoother)
    d = ndims(grids)
    if grid_ptr == sum(size(grids))-1
        grids[grid_ptr].x .= zeros(grids[grid_ptr].x)
        sol = Matrix{Vector{eltype(grids[1].x)}}(size(grids))
    else
        for idx in NotSoSimpleMultigrid.grids_at_level(size(grids),grid_ptr+1)
            child_iter = Base.Iterators.filter(i->all(i[2].>=1),enumerate([idx.-NotSoSimpleMultigrid.δ(i,d) for i in 1:d]))
            grids[idx...].b .= mean(map(i->grids[last(i)...].R[first(i)]*grids[last(i)...].b,child_iter))
        end
        sol = FMG(grids,ν₀,ν₁,ν₂,grid_ptr+1,smoother)
        for idx in NotSoSimpleMultigrid.grids_at_level(size(grids),grid_ptr)
            parent_iter = Base.Iterators.filter(i->all(i[2].<=size(grids)),enumerate([idx.+δ(i,NotSoSimpleMultigrid.d) for i in 1:d]))
            # matrix-dependent prolongation
            λ = map(i->grids[idx...].A*NotSoSimpleMultigrid.high_freq_mode(first(i),grids[idx...].sz),parent_iter)
            λ² = broadcast(i->broadcast(j->j^2,i),λ)
            ω = map(i->λ²[i]./sum(λ²),1:length(λ)) # weight factors from [Naik, Van Rosendale]
            # FMG with cubic interpolation
            ip = map(i->NotSoSimpleMultigrid.P̃(first(i),SimpleMultigrid.Cubic(),grids[last(i)...].sz...)*grids[last(i)...].x,parent_iter) # cubic
            grids[idx...].x .= sum(map(i->ω[i].*ip[i],1:length(ω)))
        end
    end
    # V-cycling
    for i in 1:ν₀
        NotSoSimpleMultigrid.μ_cycle!(grids,1,ν₁,ν₂,grid_ptr,smoother)
    end
    # safety
    for idx in NotSoSimpleMultigrid.grids_at_level(size(grids),grid_ptr)
        if SimpleMultigrid.norm_of_residu(grids[idx...]) >= 1/prod(grids[idx...].sz)
            grids[idx...].x .= grids[idx...].A\grids[idx...].b # exact solve
        end
        sol[(size(grids).-idx.+one.(idx))...] = copy(grids[idx...].x) 
    end
    return sol
end
