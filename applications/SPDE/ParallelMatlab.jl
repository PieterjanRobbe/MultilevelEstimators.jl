


# load parallel computing module
using Distributed
using MATLAB
using DelimitedFiles
# add procs
addprocs(3)

@everywhere using MATLAB
@everywhere using Distributed
a=3
# function to run in parallel
@everywhere say_hi() = eval_string(string("disp('MATLAB says hi from processor ", myid(), "')"))
for id=1:10
@time pmap(i->say_hi(), workers())

end
