# solver.jl : solver for geometric brownian motion with European callback option

function geometric_brownian_motion(level::Level, ξ::Vector{T} where {T<:Real}, sampler::Sampler)
	L = level[1]

	T = sampler.data.T # end time
	r = sampler.data.r # diffusion
	σ = sampler.data.σ # standard deviation
	M = sampler.data.M # time step refinement factor
	x₀ = sampler.data.x₀ # initial condition

	nf = M^L
	hf = T/nf
	if L == 0
		for n = 1:nf
			dWf = √hf*ξ[n]
			Xf += r*Xf*hf + σ*Xf*dWf
		end
		exp(-r*T)*Pf = max(0,Xf-1)
		Pc = 0.
	else
		nc = 4^(L-1)
		hc = T/nc
		for n = 1:nc
			dWc = 0
			for m = 1:M
				dWf = √hf*ξ[(n-1)*M+m]
				Xf += r*Xf*hf + σ*Xf*dWf
			end
			Xc += r*Xc*hc + σ*Xc*dWc
		end
		Pf = exp(-r*T)*max(0,Xf-1)
		Pc = exp(-r*T)*max(0,Xc-1)
	end
	return [Pc; Pf]
end
