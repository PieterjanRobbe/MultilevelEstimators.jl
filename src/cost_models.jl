## cost_models.jl : some default cost models used in MultilevelEstimators.jl

# multilevel cost
function ml_cost(cost_model::Function,index)
    return sum([cost_model(index) for index in union(index,keys(diff(index)))])
end

# geometric cost model
function geometric_cost_model(M::N where {N<:Real},γ::T where {T<:Real},index)
    return ml_cost((i)->sum((M.*(i.+1)).^γ),index)
end
