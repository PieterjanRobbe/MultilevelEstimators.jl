# PseudoRNG is a random number generator that generates pseudo-random numbers
# in s dimensions but is of type QMCgenerator{s}, hence it is easy to apply a RandWrapper
mutable struct PseudoRNG{s} <: QMC.QMCgenerator{s}
end

function PseudoRNG(s)
	s <= 0 && throw(BoundsError("s must be positive"))
	return PseudoRNG{s}()
end

function getPoint{s,N<:Integer}(prng::PseudoRNG{s}, k::N)
	twister = MersenneTwister(10)
	return vec(rand(twister,s,k)[:,end])
end

