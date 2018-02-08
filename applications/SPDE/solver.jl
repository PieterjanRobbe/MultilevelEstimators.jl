# solver.jl : solver for 2d lognormal diffusion problem

function lognormal_diffusion(index::Index{2}, ξ::Vector{T} where {T<:Real}, sampler::Sampler)

end

lognormal_diffusion(level::Level, ξ::Vector{T} where {T<:Real}, sampler::Sampler) = lognormal_diffusion(Index(level[1],level[1]), ξ, sampler)
