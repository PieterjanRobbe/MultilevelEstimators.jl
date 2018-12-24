## FMG.jl : custom FMG method for elliptic problems
#
# Custom FMG method for elliptic problems.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

function FMG_solve(f::Function, sz::Dims, index::Index, solver::S, ::R) where {S<:AbstractSolver, R<:AbstractReuse, B<:Bool}
    if S <: MGSolver
        mg = SimpleMultigrid.MultigridMethod(f, sz, solver.cycle)
    else
        mg = NotSoSimpleMultigrid.MultigridMethod(f, sz, solver.cycle)
    end
    mg.grids[1].b .= fill(1, size(mg.grids[1].A, 1))
    sol, ν₀s = FMG!(mg, 1)
    Rn = R <: Reuse ? CartesianIndices(index.+one(index)) : 1
    view(sol, Rn), view(size.(mg.grids), Rn), ν₀s
end

function FMG!(mg::MultigridIterable{C, G}, grid_ptr::Integer) where {C, G<:AbstractVector}
    grids = mg.grids
    if grid_ptr == length(grids)
        fill!(grids[grid_ptr].x, zero(eltype(grids[grid_ptr].x))) 
        sol = Vector{Vector{eltype(grids[1].x)}}(undef, length(grids))
        ν₀s = Vector{typeof(grid_ptr)}(undef, length(grids))
    else
        grids[grid_ptr+1].b .= grids[grid_ptr].R*grids[grid_ptr].b
        sol, ν₀s = FMG!(mg, grid_ptr+1)
        grids[grid_ptr].x .= SimpleMultigrid.P(SimpleMultigrid.Cubic(), grids[grid_ptr+1].sz...)*grids[grid_ptr+1].x
    end
    ν₀ = 0
    while !converged(grids, grid_ptr) && ν₀ < 20
        SimpleMultigrid.μ_cycle!(grids, 2, mg.cycle_type.ν₁, mg.cycle_type.ν₂, grid_ptr, mg.smoother)
        ν₀ += 1
    end
    if !converged(grids, grid_ptr) # safety
        grids[grid_ptr].x .= grids[grid_ptr].A\grids[grid_ptr].b # exact solve
    end
    sol[grid_ptr] = copy(grids[grid_ptr].x)
    ν₀s[grid_ptr] = ν₀

    return sol, ν₀s
end

converged(grid::SimpleMultigrid.Grid) = SimpleMultigrid.norm_of_residu(grid) < 1/prod(size(grid))
converged(grids::Vector{<:SimpleMultigrid.Grid}, grid_ptr) = converged(grids[grid_ptr])

function FMG!(mg::MultigridIterable{C, G}, grid_ptr::Int) where {C, G<:AbstractMatrix}
    grids = mg.grids
    R = CartesianIndices(size(grids))
    I1, Iend = first(R), last(R)
    if grid_ptr == sum(Tuple(Iend-I1)) + 1
        fill!(grids[grid_ptr].x, zero(eltype(grids[grid_ptr].x)))
        sol = Matrix{Vector{eltype(grids[1].x)}}(undef, size(grids)...)
        ν₀s = Matrix{typeof(grid_ptr)}(undef, size(grids)...)
    else
        for I in NotSoSimpleMultigrid.grids_at_level(R, grid_ptr+1)
            R_child = NotSoSimpleMultigrid.child_iter(R, I1, I)
            grids[I].b .= mean(map(i->grids[last(i)].R[first(i)]*grids[last(i)].b, R_child))
        end
        sol, ν₀s = FMG!(mg, grid_ptr+1)
        for I in NotSoSimpleMultigrid.grids_at_level(R, grid_ptr)
            R_parent = NotSoSimpleMultigrid.parent_iter(R, I1, I)
            ip = map(i->NotSoSimpleMultigrid.P̃(first(i), SimpleMultigrid.Cubic(), grids[last(i)].sz...) * grids[last(i)].x, R_parent)
            c = sum(ip)
            d = SimpleMultigrid.residu(grids[I])
            α = c'*d/(c'*grids[I].A*c)
            α = isnan(α) ? one(eltype(c)) : min(1.1, max(0.7, α))
            #grids[I].x .= α * c
            grids[I].x .= c
        end
    end
        ν₀ = 0
        while !converged(grids, grid_ptr, R) && ν₀ < 20
            NotSoSimpleMultigrid.μ_cycle!(grids, 2, mg.cycle_type.ν₁, mg.cycle_type.ν₂, 1, mg.smoother)
            ν₀ += 1
        end
    for I in NotSoSimpleMultigrid.grids_at_level(R, grid_ptr)
        if !converged(grids[I]) # safety
            grids[I].x .= grids[I].A\grids[I].b # exact solve
        end
        sol[I] = copy(grids[I].x)
        ν₀s[I] = ν₀
    end

    return sol, ν₀s
end

converged(grids::Matrix{<:SimpleMultigrid.Grid}, grid_ptr, R) = all([converged(grids[I]) for I in NotSoSimpleMultigrid.grids_at_level(R, grid_ptr)]) 

# function that returns residual norm after 50 multigrid V-cycles (used to analyze performance)
function V_cycle_solve(f::Function, sz::Dims, solver::S) where {S<:AbstractSolver}
    if S <: MGSolver
        mg = SimpleMultigrid.MultigridMethod(f, sz, solver.cycle, smoother=RedBlackGaussSeidel())
    else
        mg = NotSoSimpleMultigrid.MultigridMethod(f, sz, solver.cycle, smoother=RedBlackGaussSeidel())
    end
    mg.grids[1].b .= fill(1., size(mg.grids[1].A, 1))
    push!(mg.resnorm, SimpleMultigrid.norm_of_residu(mg.grids[1]))
    for i in 1:50
        SimpleMultigrid.cycle!(mg)
        push!(mg.resnorm, SimpleMultigrid.norm_of_residu(mg.grids[1]))
    end
    mg.resnorm
end
