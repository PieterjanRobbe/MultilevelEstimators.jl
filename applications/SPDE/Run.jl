#__precompile__()

#module Run
push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE")))
push!(LOAD_PATH,(joinpath("/","home","philippeb",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE")))
using Distributed
using MultilevelEstimators
using Coupling
using Pkg
using MATLAB
#using SharedArrays
#using DistributedArrays
#push!(LOAD_PATH,(joinpath("/","home","philippe","JuliaRuns","MultilevelEstimators","applications","SPDE")))
#push!(LOAD_PATH,Pkg.dir(joinpath("MultilevelEstimators","applications","SPDE")))
#push!(LOAD_PATH,(joinpath("/","home","philippe","JuliaRuns","MultilevelEstimators","src")))
#push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","v0.6","MultilevelEstimators","applications","SPDE")))

#./julia /home/philippeb/.julia/v0.6/MultilevelEstimators/applications/SPDE/Run.jl 2>&1 | tee /home/philippeb/.julia/v0.6/MultilevelEstimators/applications/SPDE/10112018_Run_Log.txt
numberOfProcs=1

addprocs(numberOfProcs)
println("Procs added")
#@everywhere using Coupling
@everywhere push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE")))
@everywhere push!(LOAD_PATH,(joinpath("/","home","philippeb",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE")))
@everywhere using Coupling
@everywhere using MATLAB
#@everywhere using Distributed
#@everywhere using SharedArrays
#@everywhere using DistributedArrays

if(numberOfProcs>3)
pathToInterm=joinpath("/","home","philippeb",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE","data")
else
pathToInterm=joinpath("/","home","philippe",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE","data")
end
folder = string(pathToInterm,"/Interm") # for report
if(isdir(folder)==false)
mkdir(folder)
else
    rm(folder,recursive=true)
    mkdir(folder)
end



##-------------------------------------------------------------
## NON LINEAR
##----------------------------------------------------
#-------------
# Higer order elements
#------------
#In test phase on 02/04/2019

#Heterogeneous
#mlmc
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_mlmc_NL_Het_Single_HigherOrder= Coupling.init_Beam(ML(),false,false,false,true,true,true,250/4000,nb_of_warm_up_samples=40,max_index_set_param=3,MatlabSampler,folder)
#estimator=init_Beam_mlmc_NL_Het_Single_HigherOrder
#history = run(estimator,2.5e-6)

#mlqmc
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_mlqmc_NL_Het_Single_HigherOrder=Coupling.init_Beam(ML(),true,false,false,true,true,true,250/4000,nb_of_warm_up_samples=2,max_index_set_param=3,nshifts=10,MatlabSampler,folder,false)
#estimator=init_Beam_mlqmc_NL_Het_Single_HigherOrder
#history = run(estimator,2.5e-6)

#Homogeneous
#mlmc
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_mlmc_NL_Hom_Single_HigherOrder= Coupling.init_Beam(ML(),false,false,false,true,false,true,250/4000,nb_of_warm_up_samples=40,max_index_set_param=3,MatlabSampler,folder)
#estimator=init_Beam_mlmc_NL_Hom_Single_HigherOrder
#history = run(estimator,2.5e-6)

#mlqmc
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_mlqmc_NL_Hom_Single_HigherOrder=Coupling.init_Beam(ML(),true,false,false,true,false,true,250/4000,nb_of_warm_up_samples=2,max_index_set_param=3,nshifts=10,MatlabSampler,folder,false)
#estimator=init_Beam_mlqmc_NL_Hom_Single_HigherOrder
#history = run(estimator,2.5e-6)

##-------------
#-------------
# Refinement
#------------
#finised and used for Korea Paper
#In test fase on 02/04/2019
#Heterogeneous


#mlqmc
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlqmc_NL_Het_Single= Coupling.init_Beam(ML(),true,false,false,true,true,false,250/4000,nb_of_warm_up_samples=2,max_index_set_param=3,nshifts=10,MatlabSampler,folder,false)
#estimator=init_Beam_mlqmc_NL_Het_Single
#history = run(estimator,2.5e-6)

#mlmc
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlmc_NL_Het_Single= Coupling.init_Beam(ML(),false,false,false,true,true,false,250/4000,max_index_set_param=3,MatlabSampler,folder)
#estimator=init_Beam_mlmc_NL_Het_Single
#history = run(estimator,2.5e-6)



#Homogeneous
#mlmc
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlmc_NL_Hom_Single= Coupling.init_Beam(ML(),false,false,false,true,false,false,250/4000,max_index_set_param=3,MatlabSampler,folder)
#estimator=init_Beam_mlmc_NL_Hom_Single
#history = run(estimator,2.5e-6)

#In test fase on 02/04/2019
#mlqmc`
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlqmc_NL_Hom_Single= Coupling.init_Beam(ML(),true,false,false,true,false,false,250/4000,nb_of_warm_up_samples=2,max_index_set_param=3,nshifts=10,MatlabSampler,folder,false)
#estimator=init_Beam_mlqmc_NL_Hom_Single
#history = run(estimator,2.5e-6)


##-------------------------------------------------------------
##-------------------------------------------------------------
##  LINEAR Testing on 06/04/2019
##----------------------------------------------------
#-------------
# Higer order elements
#------------
#Heterogeneous
##mlmc
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_mlmc_L_Het_Single_HigherOrder= Coupling.init_Beam(ML(),false,false,false,false,true,true,250/4000,nb_of_warm_up_samples=40,max_index_set_param=4,MatlabSampler,folder)
#estimator=init_Beam_mlmc_L_Het_Single_HigherOrder
#history = run(estimator,8.0e-5)

#mlqmc
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_mlqmc_L_Het_Single_HigherOrder=Coupling.init_Beam(ML(),true,false,false,false,true,true,250/4000,nb_of_warm_up_samples=2,max_index_set_param=4,nshifts=10,MatlabSampler,folder,false)
#estimator=init_Beam_mlqmc_L_Het_Single_HigherOrder
#history = run(estimator,8.0e-5)

#Homogeneous
#mlmc
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_mlmc_L_Hom_Single_HigherOrder= Coupling.init_Beam(ML(),false,false,false,false,false,true,250/4000,nb_of_warm_up_samples=40,max_index_set_param=5,MatlabSampler,folder,false)
#estimator=init_Beam_mlmc_L_Hom_Single_HigherOrder
#history = run(estimator,3.861e-05)

#mlqmc
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_mlqmc_L_Hom_Single_HigherOrder=Coupling.init_Beam(ML(),true,false,false,false,false,true,250/4000,nb_of_warm_up_samples=2,max_index_set_param=5,nshifts=10,MatlabSampler,folder,false)
#estimator=init_Beam_mlqmc_L_Hom_Single_HigherOrder
#history = run(estimator,3.861e-05)

#-------------
# Refinement
#------------
#Heterogeneous
#mlmc
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlmc_L_Het_Single= Coupling.init_Beam(ML(),false,false,false,false,true,false,250/4000,max_index_set_param=5,MatlabSampler,folder,false)
#estimator=init_Beam_mlmc_L_Het_Single
#history = run(estimator,3.861e-05)


#mlqmc
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlqmc_L_Het_Single= Coupling.init_Beam(ML(),true,false,false,false,true,false,250/4000,nb_of_warm_up_samples=2,max_index_set_param=5,nshifts=10,MatlabSampler,folder,false)
#estimator=init_Beam_mlqmc_L_Het_Single
#history = run(estimator,3.861e-05)

#Homogeneous
#mlmc
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlmc_L_Hom_Single= Coupling.init_Beam(ML(),false,false,false,false,false,false,250/4000,max_index_set_param=5,MatlabSampler,folder,false)
#estimator=init_Beam_mlmc_L_Hom_Single
#history = run(estimator,3.861e-05)

#mlqmc
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_mlqmc_L_Hom_Single= Coupling.init_Beam(ML(),true,false,false,false,false,false,250/4000,nb_of_warm_up_samples=2,max_index_set_param=5,nshifts=10,MatlabSampler,folder,false)
#estimator=init_Beam_mlqmc_L_Hom_Single
#history = run(estimator,3.861e-05)

#############################################################
################ MONTE CARLO
######################################################
#LINEAR

#HETEROGENEOUS

#HIGH
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_MC_L_Het_Single_high=Coupling.init_Beam(SL(),false,false,false,false,true,true,250/16000,startlevel=3,MatlabSampler,folder,false)
#estimator=init_Beam_MC_L_Het_Single_high
#history = run(estimator,3.861e-05)
#REF
@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
init_Beam_MC_L_Het_Single=Coupling.init_Beam(SL(),false,false,false,false,true,false,250/4000,startlevel=3,numberoftol=10,MatlabSampler,folder,false)
estimator=init_Beam_MC_L_Het_Single
history = run(estimator,3.861e-05)



############
#HOMOGENEOUS

#HIGH
@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
init_Beam_MC_L_Hom_Single_high=Coupling.init_Beam(SL(),false,false,false,false,false,true,250/4000,startlevel=3,numberoftol=10,MatlabSampler,folder,false)
estimator=init_Beam_MC_L_Hom_Single_high
history = run(estimator,3.861e-05)
#REF
@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_READ_IN(", myid(), ")"))
init_Beam_MC_L_Hom_Single=Coupling.init_Beam(SL(),false,false,false,false,false,false,250/4000,startlevel=3,numberoftol=6,MatlabSampler,folder,false)
estimator=init_Beam_MC_L_Hom_Single
history = run(estimator,0.000110274021000000)
#######

#NONLINEAR

#HETEROGENEOUS

#HIGH
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_MC_NL_Het_Single_high=Coupling.init_Beam(SL(),false,false,false,true,true,true,250/4000,startlevel=3,numberoftol=8,MatlabSampler,folder)
#estimator=init_Beam_MC_NL_Het_Single_high
#history = run(estimator,4.22500000000000e-06)
#REF
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_MC_NL_Het_Single=Coupling.init_Beam(SL(),false,false,false,true,true,false,250/4000,startlevel=3,numberoftol=7,MatlabSampler,folder)
#estimator=init_Beam_MC_NL_Het_Single
#history = run(estimator,5.49250000000000e-06)
############
#HOMOGENEOUS

#HIGH
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_High_READ_IN(", myid(), ")"))
#init_Beam_MC_NL_Hom_Single_high=Coupling.init_Beam(SL(),false,false,false,true,false,true,250/4000,startlevel=3,numberoftol=6,MatlabSampler,folder)
#estimator=init_Beam_MC_NL_Hom_Single_high
#history = run(estimator,7.14025000000000e-06)
#REF
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_READ_IN(", myid(), ")"))
#init_Beam_MC_NL_Hom_Single=Coupling.init_Beam(SL(),false,false,false,true,false,false,250/4000,startlevel=3,numberoftol=3,MatlabSampler,folder)
#estimator=init_Beam_MC_NL_Hom_Single
#history = run(estimator,1.56871292500000e-05)
#######


##################3
######## LINEAR GAUSS POINTS
#-------------
# Higer order elements
#------------
#Heterogeneous
##mlmc
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_High_READ_IN_GP(", myid(), ")"))
#init_Beam_mlmc_L_Het_Single_HigherOrder= Coupling.init_Beam(ML(),false,false,false,false,true,true,250/4000,nb_of_warm_up_samples=40,max_index_set_param=5,MatlabSampler,folder,true)
#estimator=init_Beam_mlmc_L_Het_Single_HigherOrder
#history = run(estimator,3.861e-05)

#mlqmc
#@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_High_READ_IN_GP(", myid(), ")"))
#init_Beam_mlqmc_L_Het_Single_HigherOrder=Coupling.init_Beam(ML(),true,false,false,false,true,true,250/4000,nb_of_warm_up_samples=2,max_index_set_param=5,nshifts=10,MatlabSampler,folder,true)
#estimator=init_Beam_mlqmc_L_Het_Single_HigherOrder
#history = run(estimator,3.861e-05)

####### MC

@everywhere MatlabSampler() = eval_string(string("Solver_L_JULIA_MATLAB_High_READ_IN_GP(", myid(), ")"))
init_Beam_MC_L_Het_Single_high=Coupling.init_Beam(SL(),false,false,false,false,true,true,250/4000,startlevel=3,MatlabSampler,folder,true)
estimator=init_Beam_MC_L_Het_Single_high
history = run(estimator,3.861e-05)

#-------------------------------------------------------------
## NON LINEAR
##----------------------------------------------------
#-------------
# Higer order elements Gausspoints
#------------
#In test phase on

#Heterogeneous
#mlmc
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_High_READ_IN_GP(", myid(), ")"))
#init_Beam_mlmc_NL_Het_Single_HigherOrder= Coupling.init_Beam(ML(),false,false,false,true,true,true,250/4000,nb_of_warm_up_samples=40,max_index_set_param=3,MatlabSampler,folder,true)
#estimator=init_Beam_mlmc_NL_Het_Single_HigherOrder
#history = run(estimator,2.5e-6)

#mlqmc
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_High_READ_IN_GP(", myid(), ")"))
#init_Beam_mlqmc_NL_Het_Single_HigherOrder=Coupling.init_Beam(ML(),true,false,false,true,true,true,250/4000,nb_of_warm_up_samples=2,max_index_set_param=3,nshifts=10,MatlabSampler,folder,true)
#estimator=init_Beam_mlqmc_NL_Het_Single_HigherOrder
#history = run(estimator,2.5e-6)

####MC

#HIGH
#@everywhere MatlabSampler() = eval_string(string("Solver_NL_JULIA_MATLAB_High_READ_IN_GP(", myid(), ")"))
#init_Beam_MC_NL_Het_Single_high=Coupling.init_Beam(SL(),false,false,false,true,true,true,250/4000,startlevel=3,numberoftol=8,MatlabSampler,folder,true)
#estimator=init_Beam_MC_NL_Het_Single_high
#history = run(estimator,4.22500000000000e-06)

#end
