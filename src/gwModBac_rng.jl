# SeedGenerator is a random number generator that generates randim UInt32 numbers that can be used as seed
mutable struct SeedGenerator
	seed
end

function getPoint{N<:Integer}(rng::SeedGenerator, k::N)
	twister = MersenneTwister(rng.seed)
	return rand(twister,UInt32,k)[end]
end

