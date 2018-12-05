## FMG.jl : custom FMG method for elliptic problems
#
# Custom FMG method for elliptic problems.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

function FMG_solve(f::Function, sz::Dims, ::S, ::R) where {S<:AbstractSolver, R<:AbstractReuse}
    if S <: MGSolver
        mg = SimpleMultigrid.MultigridMethod(f, sz, V(6, 4))
    else
        mg = NotSoSimpleMultigrid.MultigridMethod(f, sz, V(6, 4))
    end
    mg.grids[1].b .= fill(1, size(mg.grids[1].A, 1))
    sol = FMG!(mg, 1)
    if R <: NoReuse
        return sol[1], size(mg.grids[1])
    else
        return sol, size.(mg.grids)
    end
end

function FMG!(mg::MultigridIterable{C, G}, grid_ptr::Int) where {C, G<:AbstractVector}
    grids = mg.grids
    if grid_ptr == length(grids)
        fill!(grids[grid_ptr].x, zero(eltype(grids[grid_ptr].x))) 
        sol = Vector{Vector{eltype(grids[1].x)}}(undef, length(grids))
    else
        grids[grid_ptr+1].b .= grids[grid_ptr].R*grids[grid_ptr].b
        sol = FMG!(mg, grid_ptr+1)
        grids[grid_ptr].x .= SimpleMultigrid.P(SimpleMultigrid.Cubic(), grids[grid_ptr+1].sz...)*grids[grid_ptr+1].x
    end
    ν₀= 0
    while SimpleMultigrid.norm_of_residu(grids[grid_ptr]) >= 1/prod(grids[grid_ptr].sz) && ν₀ < 15
        SimpleMultigrid.μ_cycle!(grids, 1, mg.cycle_type.ν₁, mg.cycle_type.ν₂, grid_ptr, mg.smoother)
        ν₀ += 1
    end
    if SimpleMultigrid.norm_of_residu(grids[grid_ptr]) >= 1/prod(grids[grid_ptr].sz) # safety
        grids[grid_ptr].x .= grids[grid_ptr].A\grids[grid_ptr].b # exact solve
    end
    sol[grid_ptr] = copy(grids[grid_ptr].x)

    return sol
end

function FMG!(mg::MultigridIterable{C, G}, grid_ptr::Int) where {C, G<:AbstractMatrix}
    grids = mg.grids
    R = CartesianIndices(size(grids))
    I1, Iend = first(R), last(R)
    if grid_ptr == sum(Tuple(Iend-I1)) + 1
        fill!(grids[grid_ptr].x, zero(eltype(grids[grid_ptr].x)))
        sol = Matrix{Vector{eltype(grids[1].x)}}(undef, size(grids)...)
    else
        for I in NotSoSimpleMultigrid.grids_at_level(R, grid_ptr+1)
            R_child = NotSoSimpleMultigrid.child_iter(R, I1, I)
            grids[I].b .= mean(map(i->grids[last(i)].R[first(i)]*grids[last(i)].b, R_child))
        end
        sol = FMG!(mg, grid_ptr+1)
        for I in NotSoSimpleMultigrid.grids_at_level(R, grid_ptr)
            R_parent = NotSoSimpleMultigrid.parent_iter(R, I1, I)
            # matrix-dependent prolongation
            λ = map(i->grids[I].A * NotSoSimpleMultigrid.high_freq_mode(first(i), grids[I].sz), R_parent)
            λ² = broadcast(i->broadcast(j->j^2, i), λ)
            ω = map(i->λ²[i]./sum(λ²),1:length(λ)) # weight factors from [Naik, Van Rosendale]
            ip = map(i->NotSoSimpleMultigrid.P̃(first(i), SimpleMultigrid.Cubic(), grids[last(i)].sz...) * grids[last(i)].x, R_parent)
            grids[I].x .= sum(map(i->ω[i].*ip[i], 1:length(ω)))
        end
    end
    ν₀= 0
    while NotSoSimpleMultigrid.norm_of_residu(grids[grid_ptr]) >= 1/prod(grids[grid_ptr].sz) && ν₀ < 15
        NotSoSimpleMultigrid.μ_cycle!(grids, 1, mg.cycle_type.ν₁, mg.cycle_type.ν₂, grid_ptr, mg.smoother)
        ν₀ += 1
    end
    if NotSoSimpleMultigrid.norm_of_residu(grids[grid_ptr]) >= 1/prod(grids[grid_ptr].sz) # safety
        grids[grid_ptr].x .= grids[grid_ptr].A\grids[grid_ptr].b # exact solve
    end
    for I in NotSoSimpleMultigrid.grids_at_level(R, grid_ptr)
        sol[I] = copy(grids[I].x)
    end

    return sol
end

# function that returns residual norm after 50 multigrid V-cycles (used to analyze performance)
function V_cycle_solve(f::Function, sz::Dims, ::S) where {S<:AbstractSolver}
    if S <: MGSolver
        mg = SimpleMultigrid.MultigridMethod(f, sz, V(6,4), damping=damping)
    else
        mg = NotSoSimpleMultigrid.MultigridMethod(f, sz, V(6,4), damping=damping)
    end
    mg.grids[1].b .= fill(1, size(mg.grids[1].A, 1))
	push!(mg.resnorm, SimpleMultigrid.norm_of_residu(mg.grids[1]))
	for i in 1:100
		SimpleMultigrid.cycle!(mg)
		push!(mg.resnorm, SimpleMultigrid.norm_of_residu(mg.grids[1]))
	end
	mg.resnorm
end
