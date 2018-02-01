## random_field_samplers : gaussian random field samplers 

## GaussianRandomFieldSampler ##
"""
GaussianRandomFieldSampler

Abstract Gaussian random field sampler type.
"""
abstract type GaussianRandomFieldSampler end

## EmptyGaussianRandomFieldSampler ## 
"""
EmptyGaussianRandomFieldSampler

Singleton type of an empty Gaussian random field sampler.
"""
struct EmptyGaussianRandomFieldSampler <: GaussianRandomFieldSampler
end
