module KLEigenvals


using   GaussianRandomFields,PyPlot

p=2
nterms=1000
nq=Int64(ceil(sqrt(3*nterms)))
exp_field = GaussianRandomFields.Exponential(0.3,σ=1,p=p)

cov = CovarianceFunction(2,exp_field)
vx = 0:0.01:1
vy= 0:0.01:1
grfs = GaussianRandomField(cov,KarhunenLoeve(nterms),vx,vy,quad=GaussLegendre(),nq=nq)



println(grfs)
println(grfs.data.eigenval.^2)



end
