push!(LOAD_PATH,".")
using BrownianBridges

n1 = 100
n2 = 90
bridge1 = create_bridge(n1)
bridge2 = create_bridge(n2)
xi = randn(max(n1,n2))
dB1 = brownian_bridgify(bridge1,xi)
dB2 = brownian_bridgify(bridge2,xi)
B1 = zeros(dB1)
B1[1] = 1.
for i = 2:n1
	dB = sqrt(1/n1)*dB1[i]
	B1[i] = B1[i-1] + 0.75*B1[i-1]/n1 + 0.2*B1[i-1]*dB
end
B2 = zeros(dB2)
B2[1] = 1.
for i = 2:n2
	dB = sqrt(1/n2)*dB2[i]
	B2[i] = B2[i-1] + 0.75*B2[i-1]/n1 + 0.2*B2[i-1]*dB
end

using PyPlot

t1 = linspace(0,1,n1)
t2 = linspace(0,1,n2)

figure("brownian bridge")
plot(t1,B1)
plot(t2,B2)
show()
