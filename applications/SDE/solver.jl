# solver.jl : solver for geometric brownian motion with different callback options

function geometric_brownian_motion(level::Level,ξ::Vector{T} where {T<:Real})
	L = level[1]

	T = 1
	r = 0.05
	σ = 0.2
	x₀ = 1

	nf = 4^L
	hf = T/nf
	if L == 0
		for n = 1:nf
			dWf = √hf*ξ[n]
			Xf += r*Xf*hf + σ*Xf*dWf
		end
		exp(-r*T)*Pf = max(0,Xf-1)
		Pc = 0.
	else
		# TODO split dWf and dWC ==> for Brownian bridge (first n components are the most important ones)
		nc = 4^(L-1)
		hc = T/nc
		for n = 1:nc
			dWc = 0
			for m = 1:M
				dWf = √hf*ξ[(n-1)*nc+m]
				Xf += r*Xf*hf + σ*Xf*dWf
			end
			Xc += r*Xc*hc + σ*Xc*dWc
		end
		Pf = exp(-r*T)*max(0,Xf-1)
		Pc = exp(-r*T)*max(0,Xc-1)
	end
	return [Pc; Pf]
end
