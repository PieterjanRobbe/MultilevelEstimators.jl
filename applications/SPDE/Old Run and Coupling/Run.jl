__precompile__()

module Run
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

addprocs(26)
println("Procs added")
#@everywhere using Coupling
@everywhere push!(LOAD_PATH,(joinpath("/","home","philippe",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE")))
@everywhere push!(LOAD_PATH,(joinpath("/","home","philippeb",".julia","packages","MultilevelEstimators","l8j9n","applications","SPDE")))
@everywhere using Coupling
@everywhere using MATLAB
#@everywhere using Distributed
#@everywhere using SharedArrays
#@everywhere using DistributedArrays

 s1=MSession()


 array = Array{MSession}(undef,0)
 push!(array,s1)


 arrayPtr = Array{Ptr{Nothing}}(undef,3)
 arrayPtr[1]=s1.ptr

 arrayPtr_end = Array{Ptr{Nothing}}(undef,3)
 arrayPtr_end[1]=s1.bufptr








#println(array)

#using Coupling
#LINEAR




#estimatorQMLMC=Coupling.init_Beam_qmlmc_L_Het_Single()
#h = run(estimatorQMLMC,5e-4)


#estimatorMLMC=Coupling.init_Beam_mlmc_L_Het_Single()
#h = run(estimatorMLMC,5e-4)

#######TEST

#estimatorMLMC=Coupling.init_Beam_qmlmc_L_Het_Single()
#h = run(estimatorMLMC,2.5e-4)

#estimatorMLMC=Coupling.init_Beam_mlqmc_L_Het_Single_test_new_fun()
#h = run(estimatorMLMC,2.5e-5)

#estimatorMLMC=Coupling.init_Beam_mlmc_L_Het_Single_test_new_fun()
#h = run(estimatorMLMC,2.5e-5)

#estimatorMLMC=Coupling.init_Beam_mlmc_L_Het_Single()
#h = run(estimatorMLMC,2.5e-4)



#NON-LINEAR
#Korea paper runs
#estimatorQMLMC=Coupling.init_Beam_mlqmc_NL_Het_Single()
#h = run(estimatorQMLMC,2.5e-6)
#estimatorMLMC=Coupling.init_Beam_mlmc_NL_Het_Single()
#h = run(estimatorMLMC,2.5e-6)

#Highe order ref Working - tested on 14/03/2019
#estimatorMLMC_High=Coupling.init_Beam_mlmc_NL_Het_Single_HigherOrder()
#h = run(estimatorMLMC_High,2.5e-6)
#estimatorMLQMC_High=Coupling.init_Beam_mlqmc_NL_Het_Single_HigherOrder()
#h = run(estimatorMLQMC_High,2.5e-6)

#Not finished
#estimatorML=Coupling.init_Beam_ml_NL_Het_Single()
#h = run(estimatorML,2.5e-6)

#LINEAR
#estimatorMLMC_Lin_Hom=Coupling.init_Beam_mlmc_L_Hom_Single()
#h = run(estimatorMLMC_Lin_Hom,2.5e-4)

##############New implementation
##init functions
#init_Beam_mc=init_Beam(SL(),false,false,false,false;kwargs...)
#Should be working
#init_Beam_mlmc_NL_Hom_Single = init_Beam(ML(),false,false,false,true,false,false;kwargs...)

#working 26/03/2019
#init_Beam_mlmc_L_Hom_Single = Coupling.init_Beam(ML(),false,false,false,false,false,false,continuate=false)

#working 26/03/2019
#init_Beam_mlqmc_L_Hom_Single = Coupling.init_Beam(ML(),true,false,false,false,false,false,continuate=false,nb_of_warm_up_samples=4)

#init_Beam_mlmc_L_Het_Single= Coupling.init_Beam(ML(),false,false,false,false,true,false) #working

#To be checked if working
#init_Beam_qmlmc_L_Het_Single= init_Beam(ML(),true,false,false,false,true,false,nb_of_warm_up_samples=1;kwargs...)

#unfinished
#init_Beam_mlmc_L_Hom_Multiple= init_Beam(ML(),false,true,false,false,false,false;kwargs...)
#init_Beam_mlmc_L_Het_Multiple= init_Beam(ML(),false,true,false,false,true,false;kwargs...)


#MC test
#init_Beam_mc_NL_Het=Coupling.init_Beam(SL(),false,false,false,true,true,false,continuate=false,startlevel=3)
#estimator=init_Beam_mc_NL_Het
#history = run(estimator,2.5e-6)
#init_Beam_mc_NL_Het_Single=Coupling.init_Beam(SL(),false,false,false,true,true,true,continuate=false,startlevel=3)
#estimator=init_Beam_mc_NL_Het_Single
#history = run(estimator,2.5e-6)


##-------------------------------------------------------------
## NON LINEAR
##----------------------------------------------------
#-------------
# Higer order elements
#------------
#In test phase on 02/04/2019

#Heterogeneous
#mlmc
#init_Beam_mlmc_NL_Het_Single_HigherOrder= Coupling.init_Beam(ML(),false,false,false,true,true,true,250/4000,nb_of_warm_up_samples=40,max_index_set_param=3)
#estimator=init_Beam_mlmc_NL_Het_Single_HigherOrder
#history = run(estimator,2.5e-6)

#mlqmc
#init_Beam_mlqmc_NL_Het_Single_HigherOrder=Coupling.init_Beam(ML(),true,false,false,true,true,true,250/4000,nb_of_warm_up_samples=4,max_index_set_param=3,nshifts=10)
#estimator=init_Beam_mlqmc_NL_Het_Single_HigherOrder
#history = run(estimator,2.5e-6)

#Homogeneous
#mlmc
#init_Beam_mlmc_NL_Hom_Single_HigherOrder= Coupling.init_Beam(ML(),false,false,false,true,false,true,250/4000,nb_of_warm_up_samples=40,max_index_set_param=3)
#estimator=init_Beam_mlmc_NL_Hom_Single_HigherOrder
#history = run(estimator,2.5e-6)

#mlqmc
#init_Beam_mlqmc_NL_Hom_Single_HigherOrder=Coupling.init_Beam(ML(),true,false,false,true,false,true,250/4000,nb_of_warm_up_samples=4,max_index_set_param=3,nshifts=10)
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
#init_Beam_mlqmc_NL_Het_Single= Coupling.init_Beam(ML(),true,false,false,true,true,false,250/4000,nb_of_warm_up_samples=4,max_index_set_param=3,nshifts=10)
#estimator=init_Beam_mlqmc_NL_Het_Single
#history = run(estimator,2.5e-6)

#mlmc
init_Beam_mlmc_NL_Het_Single= Coupling.init_Beam(ML(),false,false,false,true,true,false,250/4000,max_index_set_param=3)
estimator=init_Beam_mlmc_NL_Het_Single
history = run(estimator,2.5e-6)



#Homogeneous
#mlmc
#init_Beam_mlmc_NL_Hom_Single= Coupling.init_Beam(ML(),false,false,false,true,false,false,250/4000,max_index_set_param=3)
#estimator=init_Beam_mlmc_NL_Hom_Single
#history = run(estimator,2.5e-6)

#In test fase on 02/04/2019
#mlqmc`
#init_Beam_mlqmc_NL_Hom_Single= Coupling.init_Beam(ML(),true,false,false,true,false,false,250/4000,nb_of_warm_up_samples=4,max_index_set_param=3,nshifts=10)
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
#mlmc
#init_Beam_mlmc_L_Het_Single_HigherOrder= Coupling.init_Beam(ML(),false,false,false,false,true,true,250/4000,nb_of_warm_up_samples=40,max_index_set_param=4)
#estimator=init_Beam_mlmc_L_Het_Single_HigherOrder
#history = run(estimator,8.0e-5)

#mlqmc
#init_Beam_mlqmc_L_Het_Single_HigherOrder=Coupling.init_Beam(ML(),true,false,false,false,true,true,250/4000,nb_of_warm_up_samples=4,max_index_set_param=4,nshifts=10)
#estimator=init_Beam_mlqmc_L_Het_Single_HigherOrder
#history = run(estimator,8.0e-5)

#Homogeneous
#mlmc
#init_Beam_mlmc_L_Hom_Single_HigherOrder= Coupling.init_Beam(ML(),false,false,false,false,false,true,250/4000,nb_of_warm_up_samples=40,max_index_set_param=4)
#estimator=init_Beam_mlmc_L_Hom_Single_HigherOrder
#history = run(estimator,8.0e-5)

#mlqmc
#init_Beam_mlqmc_L_Hom_Single_HigherOrder=Coupling.init_Beam(ML(),true,false,false,false,false,true,250/4000,nb_of_warm_up_samples=4,max_index_set_param=4,nshifts=10)
#estimator=init_Beam_mlqmc_L_Hom_Single_HigherOrder
#history = run(estimator,8.0e-5)

#-------------
# Refinement
#------------
#Heterogeneous
#mlmc
#init_Beam_mlmc_L_Het_Single= Coupling.init_Beam(ML(),false,false,false,false,true,false,250/4000,max_index_set_param=4,array,arrayPtr,arrayPtr_end)
#estimator=init_Beam_mlmc_L_Het_Single
#history = run(estimator,8.0e-5)

#init_Beam_mlmc_L_Het_Single= Coupling.init_Beam(ML(),false,false,false,false,true,false,250/4000,max_index_set_param=4)
#estimator=init_Beam_mlmc_L_Het_Single
#history = run(estimator,8.0e-4)

#mlqmc
#init_Beam_mlqmc_L_Het_Single= Coupling.init_Beam(ML(),true,false,false,false,true,false,250/4000,nb_of_warm_up_samples=4,max_index_set_param=4,nshifts=10)
#estimator=init_Beam_mlqmc_L_Het_Single
#history = run(estimator,8.0e-5)

#Homogeneous
#mlmc
#init_Beam_mlmc_L_Hom_Single= Coupling.init_Beam(ML(),false,false,false,false,false,false,250/4000,max_index_set_param=4)
#estimator=init_Beam_mlmc_L_Hom_Single
#history = run(estimator,8.0e-5)

#mlqmc`
#init_Beam_mlqmc_L_Hom_Single= Coupling.init_Beam(ML(),true,false,false,false,false,false,250/4000,nb_of_warm_up_samples=4,max_index_set_param=4,nshifts=10)
#estimator=init_Beam_mlqmc_L_Hom_Single
#history = run(estimator,8.0e-5)

#############################################################
################ MONTE CARLO
######################################################
#LINEAR

#HETEROGENEOUS

#HIGH
#init_Beam_MC_L_Het_Single_high=Coupling.init_Beam(SL(),false,false,false,false,true,true,250/4000,startlevel=3)
#estimator=init_Beam_MC_L_Het_Single_high
#history = run(estimator,8.00000000000000e-05)
#REF
#init_Beam_MC_L_Het_Single=Coupling.init_Beam(SL(),false,false,false,false,true,false,250/4000,startlevel=3,numberoftol=10)
#estimator=init_Beam_MC_L_Het_Single
#history = run(estimator,8.00000000000000e-05)
############
#HOMOGENEOUS

#HIGH
#init_Beam_MC_L_Hom_Single_high=Coupling.init_Beam(SL(),false,false,false,false,false,true,250/4000,startlevel=3,numberoftol=3)
#estimator=init_Beam_MC_L_Hom_Single_high
#history = run(estimator,0.000501988136000000)
#REF
#init_Beam_MC_L_Hom_Single=Coupling.init_Beam(SL(),false,false,false,false,false,false,250/4000,startlevel=3,numberoftol=3)
#estimator=init_Beam_MC_L_Hom_Single
#history = run(estimator,0.000501988136000000)
#######

#NONLINEAR

#HETEROGENEOUS

#HIGH
#init_Beam_MC_NL_Het_Single_high=Coupling.init_Beam(SL(),false,false,false,true,true,true,250/4000,startlevel=3)
#estimator=init_Beam_MC_NL_Het_Single_high
#history = run(estimator,2.0393e-5)
#REF
#init_Beam_MC_NL_Het_Single=Coupling.init_Beam(SL(),false,false,false,true,true,false,250/4000,startlevel=3)
#estimator=init_Beam_MC_NL_Het_Single
#history = run(estimator,2.0393e-5)
############
#HOMOGENEOUS

#HIGH
#init_Beam_MC_NL_Hom_Single_high=Coupling.init_Beam(SL(),false,false,false,true,false,true,250/4000,startlevel=3)
#estimator=init_Beam_MC_NL_Hom_Single_high
#history = run(estimator,2.0393e-5)
#REF
#init_Beam_MC_NL_Hom_Single=Coupling.init_Beam(SL(),false,false,false,true,false,false,250/4000,startlevel=3)
#estimator=init_Beam_MC_NL_Hom_Single
#history = run(estimator,2.0393e-5)
#######

end
