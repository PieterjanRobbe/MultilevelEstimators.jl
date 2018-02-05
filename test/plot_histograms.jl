# plot_histograms.jl : use PyPlot to generate some histograms

using MultilevelEstimators, PyPlot
const MLE = MultilevelEstimators

N = 1_000_000

## Uniform
u = UniformMCGenerator(-3*ones(1),1*ones(1))
X = [MLE.get_point(u,i)[1] for i in 1:N]
PyPlot.plt[:hist](X,50,density=true)
title("U(-3,1)")
show()

## Gaussian
n = NormalMCGenerator(-3*ones(1),1/4*ones(1))
X = [MLE.get_point(n,i)[1] for i in 1:N]
PyPlot.plt[:hist](X,50,density=true)
title("N(-3,1/4)")
show()

## TruncatedGaussian
t = TruncatedNormalMCGenerator(zeros(1),ones(1),-2*ones(1),2*ones(1))
X = [MLE.get_point(t,i)[1] for i in 1:N]
PyPlot.plt[:hist](X,50,density=true)
title("T(0,1,-2,2)")
show()

t = TruncatedNormalMCGenerator(8*ones(1),ones(1),7*ones(1),12*ones(1))
X = [MLE.get_point(t,i)[1] for i in 1:N]
PyPlot.plt[:hist](X,50,density=true)
title("T(0,1,-2,2)")
show()
