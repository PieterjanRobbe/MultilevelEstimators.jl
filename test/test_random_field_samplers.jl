## test_random_field_samplers : test for random field samplers 

## EmptyGaussianRandomFieldSampler ## 
verbose && print("testing empty gaussian random field sampler...")

egrfs = EmptyGaussianRandomFieldSampler()
@test egrfs === EmptyGaussianRandomFieldSampler()

verbose && println("done")
