using MultilevelEstimators, PyPlot


pdri=joinpath("/","home","philippe",".julia","packages","MultilevelEstimators","l8j9n","src","generating_vectors")
file=joinpath(pdri,"K_3600_32.txt")
nbofpoints=50
dims=1
v=vec(zeros(dims,1).+2)

lat=LatticeRule32(file, dims, nbofpoints)
lat_shift = ShiftedLatticeRule(lat)
lat_shift_2 = ShiftedLatticeRule(lat)

x=0:1/nbofpoints:1-1/nbofpoints
vect=zeros(dims,nbofpoints)
vect_shift=zeros(dims,nbofpoints)
vect_shift_2=zeros(dims,nbofpoints)

for id=1:nbofpoints
vect[:,id]=get_point(lat, id)
vect_shift[:,id]=get_point(lat_shift, id)
vect_shift_2[:,id]=get_point(lat_shift_2, id)

end
println(vect)
figure()
scatter(x,vect[1,:])
scatter(x,vect_shift[1,:])
scatter(x,vect_shift_2[1,:])

#figure()
#scatter(x,vect[2,:])
#figure()
#scatter(x,vect[101,:])
println(vect[1,:])
println(vect_shift[1,:])
println(vect_shift_2[1,:])
