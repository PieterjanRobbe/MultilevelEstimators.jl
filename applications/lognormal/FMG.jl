## FMG.jl : custom FMG method for elliptic problems
#
# Custom FMG method for elliptic problems.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

function FMG_solve(f::Function, sz::Dims, damping::Real)
    mg = MultigridMethod(f, sz, V(4,3), damping=damping)
    mg.grids[1].b .= fill(1, size(mg.grids[1].A, 1))
    FMG!(mg, 1)
    mg.grids[1].x
end

function FMG!(mg::MultigridIterable, grid_ptr::Int)
    grids = mg.grids
    if grid_ptr == length(grids)
        fill!(grids[grid_ptr].x, zero(eltype(grids[grid_ptr].x))) 
        sol = Vector{Vector{eltype(grids[1].x)}}(undef, length(grids))
    else
        grids[grid_ptr+1].b .= grids[grid_ptr].R*grids[grid_ptr].b
        sol = FMG!(mg, grid_ptr+1)
        grids[grid_ptr].x .= P(SimpleMultigrid.Cubic(), grids[grid_ptr+1].sz...)*grids[grid_ptr+1].x
    end
    ν₀= 0
    while norm_of_residu(grids[grid_ptr]) >= 1/prod(grids[grid_ptr].sz) && ν₀ < 15
        μ_cycle!(grids, 1, mg.cycle_type.ν₁, mg.cycle_type.ν₂, grid_ptr, mg.smoother, mg.damping)
        ν₀ += 1
    end
    if norm_of_residu(grids[grid_ptr]) >= 1/prod(grids[grid_ptr].sz) # safety
        grids[grid_ptr].x .= grids[grid_ptr].A\grids[grid_ptr].b # exact solve
    end
    sol[length(grids)-grid_ptr+1] = copy(grids[grid_ptr].x)

    return sol
end
