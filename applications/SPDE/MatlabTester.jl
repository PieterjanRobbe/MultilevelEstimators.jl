module MatlabTester

using   Interact,MATLAB,Pkg

buffer=4096


println("Matlab Tester")




Lx=2.5
Ly=0.25
he_1=0.0625
Z= 30e3
println(Lx/he_1)
println(Ly/he_1)

#   E = 200E3;

Lx=Lx*1000
Ly=Ly*1000
he=he_1*1000
nu=0.15
t=1000.0
# solve system
##TODO
#   u = mxcall(:Solver_NL_JULIA_MATLAB,1,Lx,Ly,he,Z,nu,fy,t)
println("Open Matlab Session")

s1=MSession(buffer)
put_variable(s1, :Lx, Lx)
put_variable(s1, :Ly, Ly)
put_variable(s1, :he, he)
put_variable(s1, :Z, Z)
put_variable(s1, :nu, nu)
put_variable(s1, :t, t)

@time eval_string(s1, "r = Solver_L_JULIA_MATLAB(Lx,Ly,he,Z,nu,t)")
@time u = get_mvariable(s1, :r)
@time u=jmatrix(u)
#eval_string(s1,"clear all")
#    eval_string("exit")
close(s1)

Lx=2.5
Ly=0.25
he_1=0.0625/2
Z= 30e3

println(Lx/he_1)
println(Ly/he_1)
#   E = 200E3;

Lx=Lx*1000
Ly=Ly*1000
he=he_1*1000
nu=0.15
t=1000.0
# solve system
##TODO
#   u = mxcall(:Solver_NL_JULIA_MATLAB,1,Lx,Ly,he,Z,nu,fy,t)
println("Open Matlab Session")
s1=MSession(buffer)
put_variable(s1, :Lx, Lx)
put_variable(s1, :Ly, Ly)
put_variable(s1, :he, he)
put_variable(s1, :Z, Z)
put_variable(s1, :nu, nu)
put_variable(s1, :t, t)

@time eval_string(s1, "r = Solver_L_JULIA_MATLAB(Lx,Ly,he,Z,nu,t)")
@time u = get_mvariable(s1, :r)
@time u=jmatrix(u)
#eval_string(s1,"clear all")
#    eval_string("exit")
close(s1)

Lx=2.5
Ly=0.25
he_1=0.0625/4
Z= 30e3

println(Lx/he_1)
println(Ly/he_1)
#   E = 200E3;

Lx=Lx*1000
Ly=Ly*1000
he=he_1*1000
nu=0.15
t=1000.0
# solve system
##TODO
#   u = mxcall(:Solver_NL_JULIA_MATLAB,1,Lx,Ly,he,Z,nu,fy,t)
println("Open Matlab Session")
s1=MSession(buffer)
put_variable(s1, :Lx, Lx)
put_variable(s1, :Ly, Ly)
put_variable(s1, :he, he)
put_variable(s1, :Z, Z)
put_variable(s1, :nu, nu)
put_variable(s1, :t, t)

@time eval_string(s1, "r = Solver_L_JULIA_MATLAB(Lx,Ly,he,Z,nu,t)")
@time u = get_mvariable(s1, :r)
@time u=jmatrix(u)
#eval_string(s1,"clear all")
#    eval_string("exit")
close(s1)

Lx=2.5
Ly=0.25
he_1=0.0625/8
Z= 30e3

println(Lx/he_1)
println(Ly/he_1)
#   E = 200E3;

Lx=Lx*1000
Ly=Ly*1000
he=he_1*1000
nu=0.15
t=1000.0
# solve system
##TODO
#   u = mxcall(:Solver_NL_JULIA_MATLAB,1,Lx,Ly,he,Z,nu,fy,t)
println("Open Matlab Session")
s1=MSession(buffer)
put_variable(s1, :Lx, Lx)
put_variable(s1, :Ly, Ly)
put_variable(s1, :he, he)
put_variable(s1, :Z, Z)
put_variable(s1, :nu, nu)
put_variable(s1, :t, t)

@time eval_string(s1, "r = Solver_L_JULIA_MATLAB(Lx,Ly,he,Z,nu,t)")
@time u = get_mvariable(s1, :r)
@time u=jmatrix(u)
#eval_string(s1,"clear all")
#    eval_string("exit")
close(s1)







end
