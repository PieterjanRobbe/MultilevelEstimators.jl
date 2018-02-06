module BrownianBridges

export BrownianBridge, create_bridge, brownian_bridgify

struct BrownianBridge{n,V,W}
	bI::V
	lI::V
	rI::V
	lW::W
	rW::W
	bS::W
end

function create_bridge(n::N) where {N<:Integer}
	bI = zeros(N,n); bI[1] = n
	lI = zeros(N,n)
	rI = zeros(N,n)
	lW = zeros(n); lW[1] = 1.
	rW = zeros(n); rW[1] = 1.
	bS = zeros(n); bS[1] = sqrt(n)

	map = zeros(N,n)
	map[n] = 1

	j = 1

	for i in 2:n
		while map[j] > 0
			j+=1
		end
		k = j
		while map[k] == 0
			k+=1
		end
		l = j+trunc(N,(k-1-j)/2)
		map[l] = i
		bI[i] = l
		lI[i] = j
		rI[i] = k
		denom = k + 1 - j
		numr1 = k - l
		lW[i] = numr1/denom
		numr2 = l + 1 - j
		rW[i] = numr2/denom
		bS[i] = sqrt(numr1*numr2/denom)
		j = k + 1 > n ? 1 : k +1
	end
	BrownianBridge{n,Vector{N},Vector{Float64}}(bI,lI,rI,lW,rW,bS)
end

function brownian_bridgify(bridge::BrownianBridge{n},ξ::Vector{T}) where {n,T<:Real}
	#@assert length(ξ) == n
	B = zeros(T,n)
	B[n] = bridge.bS[1]*ξ[1]

	for i in 2:n
		j = bridge.lI[i]
		k = bridge.rI[i]
		l = bridge.bI[i]
		B[i] = bridge.rW[i]*B[k] + bridge.bS[i]*ξ[i]
		B[i] += j == 1 ? 0 : bridge.lW[i]*B[j-1]
	end
	diff(prepend!(B,0.))
end

end
