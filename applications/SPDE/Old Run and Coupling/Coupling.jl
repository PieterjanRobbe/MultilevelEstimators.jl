__precompile__()

 module Coupling

 #sample_function = isa(index_set,SL) ? isHigerOrderRefinement ? (index,ξ)->Beam_NL_Het_MC_High(index, ξ,Lx,Ly, grfs[index],startlevel) : (index,ξ)->Beam_NL_Het_MC_Ref(index, ξ,Lx,Ly, grfs[index]) : isNL ? isField ? isHigerOrderRefinement ? (index,ξ) -> Beam_NL_Het_High(index, ξ,Lx,Ly, grfs[index]) :  (index,ξ) -> Beam_NL_Het(index, ξ,Lx,Ly, grfs[index]) :  isHigerOrderRefinement ? (index, ξ) -> Beam_NL_Hom_High(index, ξ,Lx,Ly, grfs[index]) : (index, ξ) -> Beam_NL_Hom(index, ξ,Lx,Ly, grfs[index])  : isField ? is_multiple_qoi ? (index,ξ,data)->Beam_L_Het_multiple(index,ξ,data,nb_of_qoi) : isHigerOrderRefinement ? (index, ξ) -> Beam_L_Het_High(index, ξ,Lx,Ly, grfs[index]) : (index, ξ) -> Beam_L_Het(index, ξ,Lx,Ly, grfs[index]) : is_multiple_qoi ? (index,ξ,data)->Beam_L_Het_multiple(index,ξ,data,nb_of_qoi) : isHigerOrderRefinement ? (index, ξ) -> Beam_L_Hom_High(index, ξ,Lx,Ly, grfs[index]) : (index, ξ) -> Beam_L_Hom(index, ξ,Lx,Ly, grfs[index])



#push!(LOAD_PATH,Pkg.dir(joinpath("Coupling","src")))
#push!(LOAD_PATH,Pkg.dir(joinpath("MultilevelEstimators","applications","SPDE")))
  using   GaussianRandomFields,Interact, FieldTransformation,MultilevelEstimators,Interpolations,Solver_L,HomogeneousNormalMatrixGen,MATLAB,Pkg
  using Distributions
  using Distributed

 # @reexport using MultilevelEstimators, GaussianRandomFields,HomogeneousNormalMatrixGen

#Distributions,
include("FieldTransformation.jl")

## import statements ##
 #import Base.getindex

##export statements
#export init_Beam_Test
## Continuation == true, samples on level 0 1 2 and then extrap
macro get_arg(key_name, default_value)
    @eval get_arg(args::Dict{Symbol, Any}, ::Val{$key_name}) = haskey(args, $key_name) ? args[$key_name] : $default_value
end

@get_arg :max_index_set_param 6

@get_arg :minpadding index->0

get_arg(args::Dict{Symbol,Any}, arg::Symbol) = get_arg(args, Val(arg))

get_arg(args::Dict{Symbol,Any}, arg::Val{T}) where T = throw(ArgumentError(string("in init_lognormal, invalid key ", T, " found")))

get_max_index_set(index_set, Int64) = get_index_set(index_set, Int64)

get_max_index_set(::SL, args) = [Level(0)]

get_max_index_set(::Union{AD, U}, args) = get_index_set(get_arg(args, :max_search_space), get_arg(args, :max_index_set_param))






function init_Beam(index_set::AbstractIndexSet, is_qmc::Bool, is_multiple_qoi::Bool, is_analyse::Bool,isNL::Bool,isField::Bool,isHigerOrderRefinement::Bool,he_start::T; corr_len::T=0.3, smoothness::T=1.0, nterms::N=101,
     max_level::N=5, nshifts::N=1,nb_of_warm_up_samples::N=40,max_index_set_param::N=10,continuate::Bool=true,startlevel::N=0,numberoftol::N=10,kwargs...) where{T<:AbstractFloat,N<:Integer,V<:MSession}
    # println(do_regression)
    # Dimensions of beam



     Lx = 2.5;
     Ly = 0.25;
     he = he_start;

     nelx=Lx/he
     nely=Ly/he

     args = Dict{Symbol,Any}(kwargs)
     args[:index_set] = index_set
     indices = get_max_index_set(index_set, max_level)
     minpadding = get_arg(args, :minpadding)
     ## Gaussian random fields ##
         coarse_dof = 1






if(isField)
                distributions = [MultilevelEstimators.Normal() for i in 1:nterms]

             p=2
                   exp_field = GaussianRandomFields.Exponential(corr_len,σ=smoothness,p=p)
                   println("P of covar equals")
                   println(p)

                   cov = CovarianceFunction(2,exp_field)


# Note fields are generated as vy x vy thus 160 x 40 example the solve transposes this


   if(isHigerOrderRefinement==false)
   # all other levels
   grfs=Dict()
   i=startlevel
   for index in indices

       j=i
       m = coarse_dof*2^i
       n = coarse_dof*2^j

       vx = 0:(he/m)/Lx:0.99999
       vy= 0:(he/n)/Ly:0.99999
      grfs[index] = GaussianRandomField(cov,KarhunenLoeve(nterms),vx,vy,quad=GaussLegendre())
      i=i+1
   end



    else
        i=0
        grfs=Dict()
        for index in indices
            j=i
            m = coarse_dof*2^i
            n = coarse_dof*2^j

            vx = 0:(he/m)/Lx:0.99999
            vy= 0:(he/n)/Ly:0.99999
            grfs[index] = GaussianRandomField(cov,KarhunenLoeve(nterms),vx,vy,quad=GaussLegendre())
        end


    end

#println(isField)
elseif(!isField)
############## for loop is working implement it everywhere
distributions = [MultilevelEstimators.Normal() for i in 1:1]
 if(isHigerOrderRefinement==false)
             grfs=Dict()
             i=startlevel
             for index in indices


                 j=i
                 m = coarse_dof*2^i
                 n = coarse_dof*2^j

                 vx = 0:(he/m)/Lx:0.99999
                 vy= 0:(he/n)/Ly:0.99999
                 grfs[index] = HomogeneousNormalMatrixGen.HomogeneousNormalMatrix(vx,vy,(vx,vy))
                 i=i+1
             end

else

    grfs=Dict()
    i=startlevel
    for index in indices


        j=0
        m = coarse_dof*2^i
        n = coarse_dof*2^j

        vx = 0:(he/m)/Lx:0.99999
        vy= 0:(he/n)/Ly:0.99999
        grfs[index] = HomogeneousNormalMatrixGen.HomogeneousNormalMatrix(vx,vy,(vx,vy))
    end


end

end
#println(grfs)
#println(grfs)


     if is_qmc
       sample_method=QMC()
     else
         sample_method=MC()
     end

      # name
    name = "Beam "
    name = is_analyse ? string(name,"analyse ") : name
    name = isa(index_set,AD) ? string(name,"A") : name
    name = isa(index_set,ML) ? string(name,"ML") : MultilevelEstimators.ndims(index_set) > 1 ? string(name,"MI") : name
    name = is_qmc ? string(name,"Q") : name
    name = string(name,"MC")
    name = isField ? string(name,"_Het") : string(name,"_Hom")
    name = isNL ? string(name,"_NonLin") : string(name,"_Lin")
    name = isHigerOrderRefinement ? string(name,"_High") : string(name,"_Ref")


    name = is_multiple_qoi ? string(name," (multiple)") : name
    nb_of_qoi = is_multiple_qoi ? Int(Lx/he*2^(max_level-1)+1) : 1
    sample_function = isa(index_set,SL) ? isNL ? isField ? isHigerOrderRefinement ? (index,ξ)->Beam_NL_Het_MC_High(index, ξ,Lx,Ly, grfs[index],startlevel) : (index,ξ)->Beam_NL_Het_MC_Ref(index, ξ,Lx,Ly, grfs[index]) : isHigerOrderRefinement ? (index,ξ)->Beam_NL_Hom_MC_High(index, ξ,Lx,Ly, grfs[index],startlevel) : (index,ξ)->Beam_NL_Hom_MC_Ref(index, ξ,Lx,Ly, grfs[index]) : isField ? isHigerOrderRefinement ? (index,ξ)->Beam_L_Het_MC_High(index, ξ,Lx,Ly, grfs[index],startlevel) : (index,ξ)->Beam_L_Het_MC_Ref(index, ξ,Lx,Ly, grfs[index]) : isHigerOrderRefinement ? (index,ξ)->Beam_L_Hom_MC_High(index, ξ,Lx,Ly, grfs[index],startlevel) : (index,ξ)->Beam_L_Hom_MC_Ref(index, ξ,Lx,Ly, grfs[index]) : isNL ? isField ? isHigerOrderRefinement ? (index,ξ) -> Beam_NL_Het_High(index, ξ,Lx,Ly, grfs[index]) :  (index,ξ) -> Beam_NL_Het(index, ξ,Lx,Ly, grfs[index]) :  isHigerOrderRefinement ? (index, ξ) -> Beam_NL_Hom_High(index, ξ,Lx,Ly, grfs[index]) : (index, ξ) -> Beam_NL_Hom(index, ξ,Lx,Ly, grfs[index])  : isField ? is_multiple_qoi ? (index,ξ,data)->Beam_L_Het_multiple(index,ξ,data,nb_of_qoi) : isHigerOrderRefinement ? (index, ξ) -> Beam_L_Het_High(index, ξ,Lx,Ly, grfs[index]) : (index, ξ) -> Beam_L_Het(index, ξ,Lx,Ly, grfs[index]) : is_multiple_qoi ? (index,ξ,data)->Beam_L_Het_multiple(index,ξ,data,nb_of_qoi) : isHigerOrderRefinement ? (index, ξ) -> Beam_L_Hom_High(index, ξ,Lx,Ly, grfs[index]) : (index, ξ) -> Beam_L_Hom(index, ξ,Lx,Ly, grfs[index])

    sample_function2 = isa(index_set,SL) ? isNL ? isField ? isHigerOrderRefinement ? (Beam_NL_Het_MC_High) : (Beam_NL_Het_MC_Ref) : isHigerOrderRefinement ? (Beam_NL_Hom_MC_High) : (Beam_NL_Hom_MC_Ref) : isField ? isHigerOrderRefinement ? (Beam_L_Het_MC_High) : (Beam_L_Het_MC_Ref) : isHigerOrderRefinement ? (Beam_L_Hom_MC_High) : (Beam_L_Hom_MC_Ref) : isNL ? isField ? isHigerOrderRefinement ? (index,ξ) -> Beam_NL_Het_High(index, ξ,Lx,Ly, grfs[index]) :  (index,ξ) -> Beam_NL_Het(index, ξ,Lx,Ly, grfs[index]) :  isHigerOrderRefinement ?  Beam_NL_Hom_High :  Beam_NL_Hom  : isField ? is_multiple_qoi ? (index,ξ,data)->Beam_L_Het_multiple(index,ξ,data,nb_of_qoi) : isHigerOrderRefinement ? (index, ξ) -> Beam_L_Het_High(index, ξ,Lx,Ly, grfs[index]) : (index, ξ) -> Beam_L_Het(index, ξ,Lx,Ly, grfs[index]) : is_multiple_qoi ? (index,ξ,data)->Beam_L_Het_multiple(index,ξ,data,nb_of_qoi) : isHigerOrderRefinement ? (index, ξ) -> Beam_L_Hom_High(index, ξ,Lx,Ly, grfs[index]) : (index, ξ) -> Beam_L_Hom(index, ξ,Lx,Ly, grfs[index])
println(sample_function)

#sample_function = (index, ξ) -> Beam_L_Het_test(index, ξ,Lx,Ly, grfs[index],MSessarray,PtrArray,PtrArray_end)

folder = string(joinpath(Pkg.dir("MultilevelEstimators"),"applications","SPDE","data",name)) # for report
if(isdir(folder)==false)
mkdir(folder)
else
    rm(folder,recursive=true)
    mkdir(folder)
end
    ## Estimator ##
    #println(folder)
if(!isa(index_set,SL))
    if(!isHigerOrderRefinement)
    γ = 2
else
    γ = 1.5
end
println("Gamma is")
println(γ)

if(is_qmc)


    MultilevelEstimators.Estimator(
    index_set, # index_set: ML, SL, TD...
    sample_method,
    sample_function,
    distributions,
    name = name, # estimator name
    folder = folder, # for report
    do_mse_splitting=false,
    nb_of_shifts=nshifts,
    nb_of_warm_up_samples=nb_of_warm_up_samples,
    max_index_set_param=max_index_set_param,
    continuate=continuate,
    cost_model = level -> 2^(γ * level[1]),
     #user_data = grfs, # GRF's
     #verbose = true, # display information
     #nb_of_qoi = nb_of_qoi, # number of qoi
     #cost_model = (index) -> geometric_cost_model(4,1.5,index), # cost model
     #sample_multiplication_factor = sample_multiplication_factor, # qmc multiplication factor
     #store_samples=false,
     )
else

    MultilevelEstimators.Estimator(
    index_set, # index_set: ML, SL, TD...
    sample_method,
    sample_function,
    distributions,
    name = name, # estimator name
    folder = folder, # for report
    do_mse_splitting=false,
    nb_of_warm_up_samples=nb_of_warm_up_samples,
    max_index_set_param=max_index_set_param,
    continuate=continuate,
    cost_model = level -> 2^(γ * level[1]),
     #user_data = grfs, # GRF's
     #verbose = true, # display information
     #nb_of_qoi = nb_of_qoi, # number of qoi
     #cost_model = (index) -> geometric_cost_model(4,1.5,index), # cost model
     #sample_multiplication_factor = sample_multiplication_factor, # qmc multiplication factor
     #store_samples=false,
     )

end
else
    MultilevelEstimators.Estimator(
    index_set, # index_set: ML, SL, TD...
    sample_method,
    sample_function,
    distributions,
    name = name, # estimator name
    folder = folder, # for report
    nb_of_warm_up_samples=nb_of_warm_up_samples,
    max_index_set_param=max_index_set_param,
    continuate=continuate,
    nb_of_tols=numberoftol
     )
 end


end

## user data ##
struct Field_Data{V}
    fields::V
end

#getindex(s::Field_Data,index::Index) = s.fields[index]

#Not fully functional yet

function Beam_NL_Het_MC_Ref(index::Index, ξ::Vector{T} where {T<:Real},Lx::Float64,Ly::Float64 , grf::GaussianRandomField)

    Error=true
    while(Error ==true)
    #    println("Start sampling")
       o(Lx,Ly,index,grf,Error)=try
    #       println("Start try catch")
        Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
        Zf=FieldTransformation.Transform_NL(Zf)
        he=Lx/size(Zf,1)
        Qf = Beam_sample_NL(Zf,Lx,Ly,he)
        dQ = Qf

        Error=false
    #    println("return no error")
        return (dQ,Qf,Error)
    catch ex
            println(ex)
            println("Error in return map caught")
            Zf=0
        Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
        Zf=FieldTransformation.Transform_NL(Zf)
        he=Lx/size(Zf,1)
        Qf = Beam_sample_NL(Zf,Lx,Ly,he)
        println("Succesfully sampled at initial level")
        dQ = Qf
        println("Start previous sample")
        Error=false
        println("Return map error ")
        return (dQ,Qf,Error)
        end
        (dQ,Qf,Error)=o(Lx,Ly,index,grf,Error)
        Error=Error
        dQ=dQ
        Qf=Qf

    end
    Error=Error
    dQ=dQ
    Qf=Qf
        return (dQ,Qf)
    end

    function Beam_NL_Het_MC_High(index::Index, ξ::Vector{T} where {T<:Real},Lx::Float64,Ly::Float64 , grf::GaussianRandomField,level::Int64)


        Error=true
        while(Error ==true)
        #    println("Start sampling")
           o(Lx,Ly,index,grf,Error,level)=try
        #       println("Start try catch")
            Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
            Zf=FieldTransformation.Transform_NL(Zf)
            he=Lx/size(Zf,1)
            Qf = Beam_sample_NL_High(Zf,Lx,Ly,he,level)
            dQ = Qf

            Error=false
        #    println("return no error")
            return (dQ,Qf,Error)
        catch ex
                println(ex)
                println("Error in return map caught")
                Zf=0
            Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
            Zf=FieldTransformation.Transform_NL(Zf)
            he=Lx/size(Zf,1)
            Qf = Beam_sample_NL_High(Zf,Lx,Ly,he,level)
            println("Succesfully sampled at initial level")
            dQ = Qf
            println("Start previous sample")
            Error=false
            println("Return map error ")
            return (dQ,Qf,Error)
            end
            (dQ,Qf,Error)=o(Lx,Ly,index,grf,Error,level)
            Error=Error
            dQ=dQ
            Qf=Qf

        end
        Error=Error
        dQ=dQ
        Qf=Qf
            return (dQ,Qf)
        end


        function Beam_NL_Hom_MC_Ref(index::Index, ξ::Vector{T} where {T<:Real},Lx::Float64,Ly::Float64 , grf::HomogeneousNormalMatrixGen.HomogeneousNormalMatrix)

            ξ_number=ξ[1]

            Error=true
            while(Error ==true)
            #    println("Start sampling")
               o(Lx,Ly,index,grf,Error)=try
            #       println("Start try catch")
                Zf = HomogeneousNormalMatrixGen.sample(grf,xi=ξ_number) # compute GRF
                Zf=FieldTransformation.Transform_NL(Zf)
                he=Lx/size(Zf,1)
                Qf = Beam_sample_NL(Zf,Lx,Ly,he)
                dQ = Qf

                Error=false
            #    println("return no error")
                return (dQ,Qf,Error)
            catch ex
                    println(ex)
                    println("Error in return map caught")
                    Zf=0
                Zf = HomogeneousNormalMatrixGen.sample(grf,xi=ξ_number) # compute GRF
                Zf=FieldTransformation.Transform_NL(Zf)
                he=Lx/size(Zf,1)
                Qf = Beam_sample_NL(Zf,Lx,Ly,he)
                println("Succesfully sampled at initial level")
                dQ = Qf
                println("Start previous sample")
                Error=false
                println("Return map error ")
                return (dQ,Qf,Error)
                end
                (dQ,Qf,Error)=o(Lx,Ly,index,grf,Error)
                Error=Error
                dQ=dQ
                Qf=Qf

            end
            Error=Error
            dQ=dQ
            Qf=Qf
                return (dQ,Qf)
            end

            function Beam_NL_Hom_MC_High(index::Index, ξ::Vector{T} where {T<:Real},Lx::Float64,Ly::Float64 , grf::HomogeneousNormalMatrixGen.HomogeneousNormalMatrix,level::Int64)

                ξ_number=ξ[1]

                Error=true
                while(Error ==true)
                #    println("Start sampling")
                   o(Lx,Ly,index,grf,Error,level)=try
                #       println("Start try catch")
                    Zf = HomogeneousNormalMatrixGen.sample(grf,xi=ξ_number) # compute GRF
                    Zf=FieldTransformation.Transform_NL(Zf)
                    he=Lx/size(Zf,1)
                    Qf = Beam_sample_NL_High(Zf,Lx,Ly,he,level)
                    dQ = Qf

                    Error=false
                #    println("return no error")
                    return (dQ,Qf,Error)
                catch ex
                        println(ex)
                        println("Error in return map caught")
                        Zf=0
                    Zf = HomogeneousNormalMatrixGen.sample(grf,xi=ξ_number) # compute GRF
                    Zf=FieldTransformation.Transform_NL(Zf)
                    he=Lx/size(Zf,1)
                    Qf = Beam_sample_NL_High(Zf,Lx,Ly,he,level)
                    println("Succesfully sampled at initial level")
                    dQ = Qf
                    println("Start previous sample")
                    Error=false
                    println("Return map error ")
                    return (dQ,Qf,Error)
                    end
                    (dQ,Qf,Error)=o(Lx,Ly,index,grf,Error,level)
                    Error=Error
                    dQ=dQ
                    Qf=Qf

                end
                Error=Error
                dQ=dQ
                Qf=Qf
                    return (dQ,Qf)
                end
#######
function Beam_L_Het_MC_Ref(index::Index, ξ::Vector{T} where {T<:Real},Lx::Float64,Ly::Float64 , grf::GaussianRandomField)


    Error=true
    while(Error ==true)
    #    println("Start sampling")
       o(Lx,Ly,index,grf,Error)=try
    #       println("Start try catch")
        Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
        Zf=FieldTransformation.Transform_L(Zf)
        he=Lx/size(Zf,1)
        Qf = Beam_sample_L(Zf,Lx,Ly,he)
        dQ = Qf

        Error=false
    #    println("return no error")
        return (dQ,Qf,Error)
    catch ex
            println(ex)
            println("Error in return map caught")
            Zf=0
        Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
        Zf=FieldTransformation.Transform_L(Zf)
        he=Lx/size(Zf,1)
        Qf = Beam_sample_L(Zf,Lx,Ly,he)
        println("Succesfully sampled at initial level")
        dQ = Qf
        println("Start previous sample")
        Error=false
        println("Return map error ")
        return (dQ,Qf,Error)
        end
        (dQ,Qf,Error)=o(Lx,Ly,index,grf,Error)
        Error=Error
        dQ=dQ
        Qf=Qf

    end
    Error=Error
    dQ=dQ
    Qf=Qf
        return (dQ,Qf)
    end

    function Beam_L_Het_MC_High(index::Index, ξ::Vector{T} where {T<:Real},Lx::Float64,Ly::Float64 , grf::GaussianRandomField,level::Int64)


        Error=true
        while(Error ==true)
        #    println("Start sampling")
           o(Lx,Ly,index,grf,Error,level)=try
        #       println("Start try catch")
            Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
            Zf=FieldTransformation.Transform_L(Zf)
            he=Lx/size(Zf,1)
            Qf = Beam_sample_L_High(Zf,Lx,Ly,he,level)
            dQ = Qf

            Error=false
        #    println("return no error")
            return (dQ,Qf,Error)
        catch ex
                println(ex)
                println("Error in return map caught")
                Zf=0
            Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
            Zf=FieldTransformation.Transform_L(Zf)
            he=Lx/size(Zf,1)
            Qf = Beam_sample_L_High(Zf,Lx,Ly,he,level)
            println("Succesfully sampled at initial level")
            dQ = Qf
            println("Start previous sample")
            Error=false
            println("Return map error ")
            return (dQ,Qf,Error)
            end
            (dQ,Qf,Error)=o(Lx,Ly,index,grf,Error,level)
            Error=Error
            dQ=dQ
            Qf=Qf

        end
        Error=Error
        dQ=dQ
        Qf=Qf
            return (dQ,Qf)
        end


        function Beam_L_Hom_MC_Ref(index::Index, ξ::Vector{T} where {T<:Real},Lx::Float64,Ly::Float64 , grf::HomogeneousNormalMatrixGen.HomogeneousNormalMatrix)

            ξ_number=ξ[1]

            Error=true
            while(Error ==true)
            #    println("Start sampling")
               o(Lx,Ly,index,grf,Error)=try
            #       println("Start try catch")
                Zf = HomogeneousNormalMatrixGen.sample(grf,xi=ξ_number) # compute GRF
                Zf=FieldTransformation.Transform_L(Zf)
                he=Lx/size(Zf,1)
                Qf = Beam_sample_L(Zf,Lx,Ly,he)
                dQ = Qf

                Error=false
            #    println("return no error")
                return (dQ,Qf,Error)
            catch ex
                    println(ex)
                    println("Error in return map caught")
                    Zf=0
                Zf = HomogeneousNormalMatrixGen.sample(grf,xi=ξ_number) # compute GRF
                Zf=FieldTransformation.Transform_L(Zf)
                he=Lx/size(Zf,1)
                Qf = Beam_sample_L(Zf,Lx,Ly,he)
                println("Succesfully sampled at initial level")
                dQ = Qf
                println("Start previous sample")
                Error=false
                println("Return map error ")
                return (dQ,Qf,Error)
                end
                (dQ,Qf,Error)=o(Lx,Ly,index,grf,Error)
                Error=Error
                dQ=dQ
                Qf=Qf

            end
            Error=Error
            dQ=dQ
            Qf=Qf
                return (dQ,Qf)
            end

            function Beam_L_Hom_MC_High(index::Index, ξ::Vector{T} where {T<:Real},Lx::Float64,Ly::Float64 , grf::HomogeneousNormalMatrixGen.HomogeneousNormalMatrix,level::Int64)

                ξ_number=ξ[1]

                Error=true
                while(Error ==true)
                #    println("Start sampling")
                   o(Lx,Ly,index,grf,Error,level)=try
                #       println("Start try catch")
                    Zf = HomogeneousNormalMatrixGen.sample(grf,xi=ξ_number) # compute GRF
                    Zf=FieldTransformation.Transform_L(Zf)
                    he=Lx/size(Zf,1)
                    Qf = Beam_sample_L_High(Zf,Lx,Ly,he,level)
                    dQ = Qf

                    Error=false
                #    println("return no error")
                    return (dQ,Qf,Error)
                catch ex
                        println(ex)
                        println("Error in return map caught")
                        Zf=0
                    Zf = HomogeneousNormalMatrixGen.sample(grf,xi=ξ_number) # compute GRF
                    Zf=FieldTransformation.Transform_L(Zf)
                    he=Lx/size(Zf,1)
                    Qf = Beam_sample_L_High(Zf,Lx,Ly,he,level)
                    println("Succesfully sampled at initial level")
                    dQ = Qf
                    println("Start previous sample")
                    Error=false
                    println("Return map error ")
                    return (dQ,Qf,Error)
                    end
                    (dQ,Qf,Error)=o(Lx,Ly,index,grf,Error,level)
                    Error=Error
                    dQ=dQ
                    Qf=Qf

                end
                Error=Error
                dQ=dQ
                Qf=Qf
                    return (dQ,Qf)
                end





#Fully functional
#02/04/2019
function Beam_NL_Het(index::Index, ξ::Vector{T} where {T<:Real},Lx::Float64,Ly::Float64 , grf::GaussianRandomField)



Error=true
while(Error ==true)
#    println("Start sampling")
   o(Lx,Ly,index,grf,Error)=try
#       println("Start try catch")
    Zf = GaussianRandomFields.sample(grf,xi=ξ) # compute GRF
    Zf=FieldTransformation.Transform_NL(Zf)
    he=Lx/size(Zf,1)
    Qf = Beam_sample_NL(Zf,Lx,Ly,he)
    dQ = Qf
    for (key,value) in diff(index)
        step = (index - key).I .+ 1

        Zc = Array(view(Zf, step[1]:step[1]:size(Zf, 1), step[1]:step[1]:size(Zf, 2)))
        he=Lx/size(Zc,1)
        Qc = Beam_sample_NL(Zc,Lx,Ly,he)
        dQ += value*Qc
    end
    Error=false
#    println("return no error")
    return (dQ,Qf,Error)
catch ex
        println(ex)
        println("Error in return map caught")
        Zf=0
    Zf = GaussianRandomFields.sample(grf,xi=ξ) # compute GRF
    Zf=FieldTransformation.Transform_NL(Zf)
    he=Lx/size(Zf,1)
    Qf = Beam_sample_NL(Zf,Lx,Ly,he)
    println("Succesfully sampled at initial level")
    dQ = Qf
    println("Start previous sample")
    for (key,value) in diff(index)
        step = (index - key).I .+ 1

        Zc = Array(view(Zf, step[1]:step[1]:size(Zf, 1), step[1]:step[1]:size(Zf, 2)))
        he=Lx/size(Zc,1)
        Qc = Beam_sample_NL(Zc,Lx,Ly,he)

        println("Succesfully sampled at previous level")
        dQ += value*Qc
    end
    Error=false
    println("Return map error ")
    return (dQ,Qf,Error)
    end
    (dQ,Qf,Error)=o(Lx,Ly,index,grf,Error)
    Error=Error
    dQ=dQ
    Qf=Qf

end
Error=Error
dQ=dQ
Qf=Qf

    #@show Qf

    # compute difference

    return (dQ,Qf)
end


function Beam_NL_Hom(index::Index, ξ::Vector{T} where {T<:Real},Lx::Float64,Ly::Float64 , grf::HomogeneousNormalMatrixGen.HomogeneousNormalMatrix)

ξ_number=ξ[1]

Error=true
while(Error ==true)
#    println("Start sampling")
   o(Lx,Ly,index,grf,Error)=try
#       println("Start try catch")
    Zf = HomogeneousNormalMatrixGen.sample(grf,ξ_number) # compute GRF
    Zf=FieldTransformation.Transform_NL(Zf)
    he=Lx/size(Zf,1)
    Qf = Beam_sample_NL(Zf,Lx,Ly,he)
    dQ = Qf
    for (key,value) in diff(index)
        step = (index - key).I .+ 1

        Zc = Array(view(Zf, step[1]:step[1]:size(Zf, 1), step[1]:step[1]:size(Zf, 2)))
        he=Lx/size(Zc,1)
        Qc = Beam_sample_NL(Zc,Lx,Ly,he)
        dQ += value*Qc
    end
    Error=false
#    println("return no error")
    return (dQ,Qf,Error)
catch ex
        println(ex)
        println("Error in return map caught")
        Zf=0
    Zf = HomogeneousNormalMatrixGen.sample(grf,ξ_number) # compute GRF
    Zf=FieldTransformation.Transform_NL(Zf)
    he=Lx/size(Zf,1)
    Qf = Beam_sample_NL(Zf,Lx,Ly,he)
    println("Succesfully sampled at initial level")
    dQ = Qf
    println("Start previous sample")
    for (key,value) in diff(index)
        step = (index - key).I .+ 1

        Zc = Array(view(Zf, step[1]:step[1]:size(Zf, 1), step[1]:step[1]:size(Zf, 2)))
        he=Lx/size(Zc,1)
        Qc = Beam_sample_NL(Zc,Lx,Ly,he)

        println("Succesfully sampled at previous level")
        dQ += value*Qc
    end
    Error=false
    println("Return map error ")
    return (dQ,Qf,Error)
    end
    (dQ,Qf,Error)=o(Lx,Ly,index,grf,Error)
    Error=Error
    dQ=dQ
    Qf=Qf

end
Error=Error
dQ=dQ
Qf=Qf

    #@show Qf

    # compute difference

    return (dQ,Qf)
end

## sample functions ##
#Working on 02/04/2019
function Beam_sample_NL(Z::Matrix{T},Lx::T,Ly::T,he::T) where {T<:Real}


 #   E = 200E3;
    nu = 0.25;
    fy = 240.0;
    t = 1.0;
    Lx=Lx*1000
    Ly=Ly*1000
    he=he*1000

    # solve system
    ##TODO
 #   u = mxcall(:Solver_NL_JULIA_MATLAB,1,Lx,Ly,he,Z,nu,fy,t)
 buffer=4096
  s1=MSession(buffer)
    put_variable(s1, :Lx, Lx)
    put_variable(s1, :Ly, Ly)
    put_variable(s1, :he, he)
    put_variable(s1, :Z, Z)
    put_variable(s1, :nu, nu)
    put_variable(s1, :fy, fy)
    put_variable(s1, :t, t)
    eval_string(s1, "r = Solver_NL_JULIA_MATLAB(Lx,Ly,he,Z,nu,fy,t);")
    u = get_mvariable(s1, :r)
    u=jmatrix(u)
#    delete(r)
    eval_string(s1,"clear all")
#    eval_string("exit")
    close(s1)
#    println("Success Matlab")
    Dispx,Dispy=Extract_matrix_Displacements_From_Vector(u[:,end],Lx,Ly,he);
#    println("Success Extrap Dir")

    Dispy=abs.(Dispy)
    Dispy=Dispy./1000
    u=0
    Dispx=0
#    clear!(:u)
#    clear!(:Dispx)
    maximalDef=findmax(Dispy)
    Dispy=0
#    clear!(:Dispy)
    GC.gc()

#   println(maximalDef[1])

     return maximalDef[1]
end

     function Beam_L_Het(index::Index, ξ::Vector{T} where {T<:Real}, Lx::Float64,Ly::Float64 , grf::GaussianRandomField)
         Error=true
         while(Error ==true)
         #    println("Start sampling")
            o(Lx,Ly,index,grf,Error)=try
         #       println("Start try catch")
             Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
             Zf=FieldTransformation.Transform_L(Zf)
             he=Lx/size(Zf,1)
             Qf = Beam_sample_L(Zf,Lx,Ly,he)
             dQ = Qf
             for (key,value) in diff(index)
                 step = (index - key).I .+ 1

                 Zc = Array(view(Zf, step[1]:step[1]:size(Zf, 1), step[1]:step[1]:size(Zf, 2)))
                 he=Lx/size(Zc,1)
                 Qc = Beam_sample_L(Zc,Lx,Ly,he)
                 dQ += value*Qc
             end
             Error=false
         #    println("return no error")
             return (dQ,Qf,Error)
         catch ex
                 println(ex)
                 println("Error in return map caught")
                 Zf=0
             Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
             Zf=FieldTransformation.Transform_L(Zf)
             he=Lx/size(Zf,1)
             Qf = Beam_sample_L(Zf,Lx,Ly,he)
             println("Succesfully sampled at initial level")
             dQ = Qf
             println("Start previous sample")
             for (key,value) in diff(index)
                 step = (index - key).I .+ 1

                 Zc = Array(view(Zf, step[1]:step[1]:size(Zf, 1), step[1]:step[1]:size(Zf, 2)))
                 he=Lx/size(Zc,1)
                 Qc = Beam_sample_L(Zc,Lx,Ly,he)

                 println("Succesfully sampled at previous level")
                 dQ += value*Qc
             end
             Error=false
             println("Return map error ")
             return (dQ,Qf,Error)
             end
             (dQ,Qf,Error)=o(Lx,Ly,index,grf,Error)
             Error=Error
             dQ=dQ
             Qf=Qf

         end
         Error=Error
         dQ=dQ
         Qf=Qf

             #@show Qf

             # compute difference

             return (dQ,Qf)
     end



     function Beam_L_Het_test(index::Index, ξ::Vector{T} where {T<:Real}, Lx::Float64,Ly::Float64 , grf::GaussianRandomField,MSessarray::AbstractArray{MSession},PtrArray::AbstractArray{Ptr{Nothing}},PtrArray_end::AbstractArray{Ptr{Nothing}})



         Error=true
         while(Error ==true)
         #    println("Start sampling")
            o(Lx,Ly,index,grf,Error)=try
         #       println("Start try catch")
             Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
             Zf=FieldTransformation.Transform_L(Zf)
             he=Lx/size(Zf,1)
             Qf = Beam_sample_L_test(Zf,Lx,Ly,he,MSessarray,PtrArray,PtrArray_end)
             dQ = Qf
             for (key,value) in diff(index)
                 step = (index - key).I .+ 1

                 Zc = Array(view(Zf, step[1]:step[1]:size(Zf, 1), step[1]:step[1]:size(Zf, 2)))
                 he=Lx/size(Zc,1)
                 Qc = Beam_sample_L_test(Zc,Lx,Ly,he,MSessarray,PtrArray,PtrArray_end)
                 dQ += value*Qc
             end
             Error=false
         #    println("return no error")
             return (dQ,Qf,Error)
         catch ex
                 println(ex)
                 println("Error in return map caught")
                 Zf=0
             Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
             Zf=FieldTransformation.Transform_L(Zf)
             he=Lx/size(Zf,1)
             Qf = Beam_sample_L_test(Zf,Lx,Ly,he,MSessarray,PtrArray,PtrArray_end)
             println("Succesfully sampled at initial level")
             dQ = Qf
             println("Start previous sample")
             for (key,value) in diff(index)
                 step = (index - key).I .+ 1

                 Zc = Array(view(Zf, step[1]:step[1]:size(Zf, 1), step[1]:step[1]:size(Zf, 2)))
                 he=Lx/size(Zc,1)
                 Qc = Beam_sample_L_test(Zc,Lx,Ly,he,MSessarray,PtrArray,PtrArray_end)

                 println("Succesfully sampled at previous level")
                 dQ += value*Qc
             end
             Error=false
             println("Return map error ")
             return (dQ,Qf,Error)
             end
             (dQ,Qf,Error)=o(Lx,Ly,index,grf,Error)
             Error=Error
             dQ=dQ
             Qf=Qf

         end
         Error=Error
         dQ=dQ
         Qf=Qf

             #@show Qf

             # compute difference

             return (dQ,Qf)
     end

     function Beam_L_Het_High(index::Index, ξ::Vector{T} where {T<:Real}, Lx::Float64,Ly::Float64 , grf::GaussianRandomField)
         Error=true
         while(Error ==true)
         #    println("Start sampling")
            o(Lx,Ly,index,grf,Error)=try
         #       println("Start try catch")
             Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
             Zf=FieldTransformation.Transform_L(Zf)
             he=Lx/size(Zf,1)
             Qf = Beam_sample_L_High(Zf,Lx,Ly,he,index[1])
             dQ = Qf
             for (key,value) in diff(index)
                 step = (index - key).I .+ 1

                 Zc = Zf
                 he=Lx/size(Zc,1)
                 Qc = Beam_sample_L_High(Zc,Lx,Ly,he,(index[1]-1))
                 dQ += value*Qc
             end
             Error=false
         #    println("return no error")
             return (dQ,Qf,Error)
         catch ex
                 println(ex)
                 println("Error in return map caught")
                 Zf=0
             Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
             Zf=FieldTransformation.Transform_L(Zf)
             he=Lx/size(Zf,1)
             Qf = Beam_sample_L_High(Zf,Lx,Ly,he,index[1])
             println("Succesfully sampled at initial level")
             dQ = Qf
             println("Start previous sample")
             for (key,value) in diff(index)
                 step = (index - key).I .+ 1

                 Zc = Zf
                 he=Lx/size(Zc,1)
                 Qc = Beam_sample_L_High(Zc,Lx,Ly,he,(index[1]-1))

                 println("Succesfully sampled at previous level")
                 dQ += value*Qc
             end
             Error=false
             println("Return map error ")
             return (dQ,Qf,Error)
             end
             (dQ,Qf,Error)=o(Lx,Ly,index,grf,Error)
             Error=Error
             dQ=dQ
             Qf=Qf

         end
         Error=Error
         dQ=dQ
         Qf=Qf

             #@show Qf

             # compute difference

             return (dQ,Qf)
     end


#To be checked
     function Beam_L_Het_multiple(index::Index, ξ::Vector{T} where {T<:Real},  data::Field_Data, nb_qoi)

                  Lx = 2.5;
                  Ly = 0.25;

     # extract grf
     grf = data[index]
#     @show index[1]
     he = 250/(2000*(2^index[1]));
     println(he)
     println(index)
     println(index[1])

     # solve


     Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF


     #Transform to Gamma Random Field
     Zf=FieldTransformation.Transform_L(Zf)

     Qf = Beam_sample_L_multiple(Zf,Lx,Ly,he,nb_qoi)

 #    @show Qf

     # compute difference
     dQ = Qf
     for (key,value) in diff(index)
#         figure()
#         surf(Zf)
         Zc = interpolate_field(data[index].pts,data[key].pts,Zf) # interpolation of fine grid GRF
#         figure()
#         surf(Zc)
         he = 250/(2000*(2^(index[1]-1)));
         Qc = Beam_sample_L_multiple(Zc,Lx,Ly,he,nb_qoi)
         dQ += value*Qc
     end

     return (dQ,Qf)
     end


         function Beam_L_Hom(index::Index, ξ::Vector{T} where {T<:Real}, Lx::Float64,Ly::Float64 ,grf::HomogeneousNormalMatrixGen.HomogeneousNormalMatrix)

        #          Lx = 2.5;
        #          Ly = 0.25;


    # he = 250/(4000*(2^.index.I));

     # solve

     ξ_number=ξ[1]



     Error=true
             while(Error ==true)
             #    println("Start sampling")
                o(Lx,Ly,index,grf,Error)=try
             #       println("Start try catch")
             Zf = HomogeneousNormalMatrixGen.sample(grf,ξ_number) # compute GRF
             Zf=FieldTransformation.Transform_L(Zf)
                 he=Lx/size(Zf,1)
                 Qf = Beam_sample_L(Zf,Lx,Ly,he)
                 dQ = Qf
                 for (key,value) in diff(index)
                     step = (index - key).I .+ 1

                     Zc = Array(view(Zf, step[1]:step[1]:size(Zf, 1), step[1]:step[1]:size(Zf, 2)))
                     he=Lx/size(Zc,1)
                     Qc = Beam_sample_L(Zc,Lx,Ly,he)
                     dQ += value*Qc
                 end
                 Error=false
             #    println("return no error")
                 return (dQ,Qf,Error)
             catch ex
                     println(ex)
                     println("Error in return map caught")
                     Zf=0
                     Zf = HomogeneousNormalMatrixGen.sample(grf,ξ_number) # compute GRF
                     Zf=FieldTransformation.Transform_L(Zf)
                 he=Lx/size(Zf,1)
                 Qf = Beam_sample_L(Zf,Lx,Ly,he)
                 println("Succesfully sampled at initial level")
                 dQ = Qf
                 println("Start previous sample")
                 for (key,value) in diff(index)
                     step = (index - key).I .+ 1

                     Zc = Array(view(Zf, step[1]:step[1]:size(Zf, 1), step[1]:step[1]:size(Zf, 2)))
                     he=Lx/size(Zc,1)
                     Qc = Beam_sample_L(Zc,Lx,Ly,he)

                     println("Succesfully sampled at previous level")
                     dQ += value*Qc
                 end
                 Error=false
                 println("Return map error ")
                 return (dQ,Qf,Error)
                 end
                 (dQ,Qf,Error)=o(Lx,Ly,index,grf,Error)
                 Error=Error
                 dQ=dQ
                 Qf=Qf

             end
             Error=Error
             dQ=dQ
             Qf=Qf

                 #@show Qf

                 # compute difference

                 return (dQ,Qf)
     end

     function Beam_L_Hom_High(index::Index, ξ::Vector{T} where {T<:Real}, Lx::Float64,Ly::Float64 , grf::HomogeneousNormalMatrixGen.HomogeneousNormalMatrix)

     #          Lx = 2.5;
     #          Ly = 0.25;


     # he = 250/(4000*(2^.index.I));

     # solve

     ξ_number=ξ[1]



     Error=true
         while(Error ==true)
         #    println("Start sampling")
            o(Lx,Ly,index,grf,Error)=try
         #       println("Start try catch")
         Zf = HomogeneousNormalMatrixGen.sample(grf,ξ_number) # compute GRF
         Zf=FieldTransformation.Transform_L(Zf)
             he=Lx/size(Zf,1)
             Qf = Beam_sample_L_High(Zf,Lx,Ly,he,index[1])
             dQ = Qf
             for (key,value) in diff(index)
                 step = (index - key).I .+ 1

                 Zc = Zf
                 he=Lx/size(Zc,1)
                 Qc = Beam_sample_L_High(Zc,Lx,Ly,he,(index[1]-1))
                 dQ += value*Qc
             end
             Error=false
         #    println("return no error")
             return (dQ,Qf,Error)
         catch ex
                 println(ex)
                 println("Error in return map caught")
                 Zf=0
                 Zf = HomogeneousNormalMatrixGen.sample(grf,ξ_number) # compute GRF
                 Zf=FieldTransformation.Transform_L(Zf)
             he=Lx/size(Zf,1)
             Qf = Beam_sample_L_High(Zf,Lx,Ly,he,index[1])
             println("Succesfully sampled at initial level")
             dQ = Qf
             println("Start previous sample")
             for (key,value) in diff(index)
                 step = (index - key).I .+ 1

                 Zc = Zf
                 he=Lx/size(Zc,1)
                 Qc = Beam_sample_L_High(Zc,Lx,Ly,he,(index[1]-1))

                 println("Succesfully sampled at previous level")
                 dQ += value*Qc
             end
             Error=false
             println("Return map error ")
             return (dQ,Qf,Error)
             end
             (dQ,Qf,Error)=o(Lx,Ly,index,grf,Error)
             Error=Error
             dQ=dQ
             Qf=Qf

         end
         Error=Error
         dQ=dQ
         Qf=Qf

             #@show Qf

             # compute difference

             return (dQ,Qf)
     end

#To be checked
     function Beam_L_Hom_multiple(index::Index, ξ::Vector{T} where {T<:Real}, data::Field_Data, nb_qoi)

                  Lx = 2.5;
                  Ly = 0.25;

     # extract grf
     grf = data[index]
#     @show index[1]
     he = 250/(2000*(2^index[1]));

     # solve


     Zf = HomogeneousNormalMatrixGen.sample(grf) # compute GRF



     #Transform to Gamma Random Field
     Zf=FieldTransformation.Transform_L(Zf)


     Qf = Beam_sample_L_multiple(Zf,Lx,Ly,he,nb_qoi)

 #    @show Qf

     # compute difference
     dQ = Qf
     for (key,value) in diff(index)
#         figure()
#         surf(Zf)
         Zc = interpolate_field(data[index].pts,data[key].pts,Zf) # interpolation of fine grid GRF
#         figure()
#         surf(Zc)
         he = 250/(2000*(2^(index[1]-1)));
         Qc = Beam_sample_L_multiple(Zc,Lx,Ly,he,nb_qoi)
         dQ += value*Qc
     end

     return (dQ,Qf)
     end

     function Beam_sample_L_multiple(Z::Matrix{T},Lx::T,Ly::T,he::T,nb_qoi::Int64) where {T<:Real}


     nu=0.15;



     # solve system
     ##TODO
#     u = mxcall(:Solver_L_JULIA_MATLAB,1,Z,he,Lx,Ly,nu)
     u=Solver_L.Interface(Z,he,Lx,Ly,nu)
     u=vec(u)
     Dispx,Dispy=Extract_matrix_Displacements_From_Vector(u,Lx,Ly,he);
     Dispy=abs.(Dispy)
     #maximalDef=findmax(Dispy)


     Dispy=Dispy[1,:]
     xvec=0:1:size(Dispy)[1]-1
     #println(size(xvec))

     #println(size(Dispy))
     #println(size(Dispy)[1]<nb_qoi)
     while size(Dispy)[1]<nb_qoi
     Dispy=FieldTransformation.interp1_equi(xvec,Dispy,2)
     xvec=0:1:size(Dispy)[1]-1
     end
     while(size(Dispy)[1]>nb_qoi)

     Dispy=FieldTransformation.DropValues(Dispy)
     end


      u=0
      Dispx=0
      clear!(:u)
      clear!(:Dispx)
      gc()

      return Dispy

     end


     #Working 01/04/2019
     ## sample functions ##
     function Beam_sample_L(Z::Matrix{T},Lx::T,Ly::T,he::T) where {T<:Real}

     nu=0.15;
     t = 1000.0;
     Lx=Lx*1000
     Ly=Ly*1000
     he=he*1000



    # buffer=4096
@time      s1=MSession()
           put_variable(s1, :Lx, Lx)
           put_variable(s1, :Ly, Ly)
           put_variable(s1, :he, he)
           put_variable(s1, :Z, Z)
           put_variable(s1, :nu, nu)
          put_variable(s1, :t, t)
          eval_string(s1, "r = Solver_L_JULIA_MATLAB(Lx,Ly,he,Z,nu,t);")
          u = get_mvariable(s1, :r)
         u=jmatrix(u)
     #    delete(r)
         eval_string(s1,"clear all")
     #    eval_string("exit")
@time         close(s1)
     #    println("Success Matlab")
         Dispx,Dispy=Extract_matrix_Displacements_From_Vector(u[:,end],Lx,Ly,he);
     #    println("Success Extrap Dir")

         Dispy=abs.(Dispy)
         Dispy=Dispy./1000
         u=0
         Dispx=0
     #    clear!(:u)
     #    clear!(:Dispx)
         maximalDef=findmax(Dispy)
         Dispy=0
     #    clear!(:Dispy)
         GC.gc()

     #   println(maximalDef[1])

          return maximalDef[1]
     end

     function Beam_sample_L_test(Z::Matrix{T},Lx::T,Ly::T,he::T,ArrayMSession::AbstractArray{MSession},PtrArray::AbstractArray{Ptr{Nothing}},PtrArray_end::AbstractArray{Ptr{Nothing}}) where {T<:Real}

     nu=0.15;
     t = 1000.0;
     Lx=Lx*1000
     Ly=Ly*1000
     he=he*1000
    #println(ArrayMSession[myid()])

     # buffer=4096
           s1=ArrayMSession[myid()]
           s1Ptr=PtrArray[myid()]
           s1Ptr_end=PtrArray_end[myid()]
           #println(myid())
           #println(s1Ptr)
           #println(s1Ptr_end)

          # println(11)
            s1.ptr=s1Ptr;
            s1.bufptr=pointer(s1.buffer);

         #  println("PointersAssigned")

           put_variable(s1, :Lx, Lx)
           put_variable(s1, :Ly, Ly)

           put_variable(s1, :he, he)
           put_variable(s1, :Z, Z)
           put_variable(s1, :nu, nu)
          put_variable(s1, :t, t)
          eval_string(s1, "r = Solver_L_JULIA_MATLAB(Lx,Ly,he,Z,nu,t);")
          u = get_mvariable(s1, :r)
         u=jmatrix(u)
     #    delete(r)
         eval_string(s1,"clear all")
     #    eval_string("exit")
#     @time         close(s1)
     #    println("Success Matlab")
         Dispx,Dispy=Extract_matrix_Displacements_From_Vector(u[:,end],Lx,Ly,he);
     #    println("Success Extrap Dir")

         Dispy=abs.(Dispy)
         Dispy=Dispy./1000
         u=0
         Dispx=0
     #    clear!(:u)
     #    clear!(:Dispx)
         maximalDef=findmax(Dispy)
         Dispy=0
     #    clear!(:Dispy)
         GC.gc()

     #   println(maximalDef[1])

          return maximalDef[1]
     end


     function Beam_sample_L_High(Z::Matrix{T},Lx::T,Ly::T,he::T,level::Int64) where {T<:Real}

     nu=0.15;
     t = 1000.0;
     Lx=Lx*1000
     Ly=Ly*1000
     he=he*1000


    buffer=4096
     s1=MSession(buffer)
         put_variable(s1, :Lx, Lx)
         put_variable(s1, :Ly, Ly)
         put_variable(s1, :he, he)
         put_variable(s1, :Z, Z)
         put_variable(s1, :nu, nu)
         put_variable(s1, :t, t)
         put_variable(s1, :level, level)
         eval_string(s1, "r = Solver_L_JULIA_MATLAB_High(Lx,Ly,he,Z,nu,t,level);")
         u = get_mvariable(s1, :r)
         u=jscalar(u)
     #    delete(r)
     eval_string(s1,"clear all")
 close(s1)
 u=u/1000;
 GC.gc()
 #   println(maximalDef[1])
  return u
     end

   #Working 03/04/2019
     function Beam_NL_Het_High(index::Index, ξ::Vector{T} where {T<:Real},Lx::Float64,Ly::Float64 , grf::GaussianRandomField)
     #println("Start Beam NL Het High")

     #println(index[1])
     Error=true
     while(Error ==true)
     #    println("Start sampling")
     o(Lx,Ly,index,grf,Error)=try
     #       println("Start try catch")
     Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
     Zf=FieldTransformation.Transform_NL(Zf)
     he=Lx/size(Zf,1)
    # println(size(Zf))

     Qf = Beam_sample_NL_High(Zf,Lx,Ly,he,index[1])
     dQ = Qf
     for (key,value) in diff(index)
        step = (index - key).I .+ 1
        Zc=Zf
    #    Zc = Array(view(Zf, step[1]:step[1]:size(Zf, 1), step[1]:step[1]:size(Zf, 2)))
    #    he=Lx/size(Zc,1)

         Qc = Beam_sample_NL_High(Zc,Lx,Ly,he,(index[1]-1))
         dQ += value*Qc
     end
     Error=false
     #    println("return no error")
     return (dQ,Qf,Error)
     catch ex
         println(ex)
         println("Error in return map caught")
         Zf=0
     Zf = GaussianRandomFields.sample(grf,xi=ξ[1:randdim(grf)]) # compute GRF
     Zf=FieldTransformation.Transform_NL(Zf)
     he=Lx/size(Zf,1)
     Qf = Beam_sample_NL_High(Zf,Lx,Ly,he,index[1])
     println("Succesfully sampled at initial level")
     dQ = Qf
     println("Start previous sample")
     for (key,value) in diff(index)
        step = (index - key).I .+ 1
        Zc=Zf
        #Zc = Array(view(Zf, step[1]:step[1]:size(Zf, 1), step[1]:step[1]:size(Zf, 2)))
        #he=Lx/size(Zc,1)
         Qc = Beam_sample_NL_High(Zc,Lx,Ly,he,(index[1]-1))
         println("Succesfully sampled at previous level")
         dQ += value*Qc
     end
     Error=false
     println("Return map error ")
     return (dQ,Qf,Error)
     end
     (dQ,Qf,Error)=o(Lx,Ly,index,grf,Error)
     Error=Error
     dQ=dQ
     Qf=Qf
     #    println(dQ)
     #    println(Qf)
     end
     Error=Error
     dQ=dQ
     Qf=Qf

     #@show Qf

     # compute difference

     return (dQ,Qf)
     end

     function Beam_NL_Hom_High(index::Index, ξ::Vector{T} where {T<:Real},Lx::Float64,Ly::Float64 , grf::HomogeneousNormalMatrixGen.HomogeneousNormalMatrix)

     #println(index[1])
     ξ_number=ξ[1]

     Error=true
     while(Error ==true)
     #    println("Start sampling")
     o(Lx,Ly,index,grf,Error)=try
     #       println("Start try catch")


     Zf = HomogeneousNormalMatrixGen.sample(grf,ξ_number) # compute GRF

     Zf=FieldTransformation.Transform_NL(Zf)
     he=Lx/size(Zf,1)
    # println(size(Zf))
     Qf = Beam_sample_NL_High(Zf,Lx,Ly,he,index[1])
     dQ = Qf
     for (key,value) in diff(index)
        step = (index - key).I .+ 1
        Zc=Zf
    #    Zc = Array(view(Zf, step[1]:step[1]:size(Zf, 1), step[1]:step[1]:size(Zf, 2)))
    #    he=Lx/size(Zc,1)
         Qc = Beam_sample_NL_High(Zc,Lx,Ly,he,(index[1]-1))
         dQ += value*Qc
     end
     Error=false
     #    println("return no error")
     return (dQ,Qf,Error)
     catch ex
         println(ex)
         println("Error in return map caught")
         Zf=0
     Zf = HomogeneousNormalMatrixGen.sample(grf,ξ_number) # compute GRF
     Zf=FieldTransformation.Transform_NL(Zf)
     he=Lx/size(Zf,1)
     Qf = Beam_sample_NL_High(Zf,Lx,Ly,he,index[1])
     println("Succesfully sampled at initial level")
     dQ = Qf
     println("Start previous sample")
     for (key,value) in diff(index)
        step = (index - key).I .+ 1
        Zc=Zf
        #Zc = Array(view(Zf, step[1]:step[1]:size(Zf, 1), step[1]:step[1]:size(Zf, 2)))
        #he=Lx/size(Zc,1)
         Qc = Beam_sample_NL_High(Zc,Lx,Ly,he,(index[1]-1))
         println("Succesfully sampled at previous level")
         dQ += value*Qc
     end
     Error=false
     println("Return map error ")
     return (dQ,Qf,Error)
     end
     (dQ,Qf,Error)=o(Lx,Ly,index,grf,Error)
     Error=Error
     dQ=dQ
     Qf=Qf
     #    println(dQ)
     #    println(Qf)
     end
     Error=Error
     dQ=dQ
     Qf=Qf

     #@show Qf

     # compute difference

     return (dQ,Qf)
     end

     ## sample functions ##
     #Updated 02/04/2019 NOT TESTED YET
     function Beam_sample_NL_High(Z::Matrix{T},Lx::T,Ly::T,he::T,level::Int64) where {T<:Real}
         nu = 0.25;
         fy = 240.0;
         t = 1.0;
         Lx=Lx*1000
         Ly=Ly*1000
         he=he*1000


     # solve system
     #Initiate Matlab link
     buffer=4096
      s1=MSession(buffer)
           put_variable(s1, :Lx, Lx)
     put_variable(s1, :Ly, Ly)
     put_variable(s1, :he, he)
     put_variable(s1, :Z, Z)
     put_variable(s1, :nu, nu)
     put_variable(s1, :fy, fy)
     put_variable(s1, :t, t)
     put_variable(s1, :level, level)
     eval_string(s1, "r = Solver_NL_JULIA_MATLAB_High(Lx,Ly,he,Z,nu,fy,t,level);")
     u = get_mvariable(s1, :r)

#     u=jmatrix(u)
     u=jscalar(u)
     #    delete(r)
     eval_string(s1,"clear all")
     close(s1)
     u=abs(u/1000);
     GC.gc()
     #   println(maximalDef[1])
      return u
     end



function Extract_matrix_Displacements_From_Vector(Points::Array{Float64,1},Lx::Float64,Ly::Float64,he::Float64)
nelx=Int64(Lx/he)
nely=Int64(Ly/he)
u=Points
#matxDir::Array{Float64,2}
#matyDir::Array{Float64,2}
#Matx::Array{Float64,2}
#Maty::Array{Float64,2}
        matxDir=zeros(nely+1,nelx+1);
        matyDir=zeros(nely+1,nelx+1);
        x=1;  y=1;  z=2;  t=1;v=1;
        while(x<=(nelx+1))
            y=1;
            while(y<=(nely+1))
                    matxDir[y,x]=u[t];
                    matyDir[y,x]=u[z];
                    z=z+2;  t=t+2;

                v=v+1;
                y=y+1;
            end
            x=x+1;
        end
        Matx=matxDir;
        Maty=matyDir;

        return Matx,Maty

end

function interpolate_field(pts_fine,pts_coarse,Z::Matrix{T}) where {T<:Real}
    itp = interpolate(pts_fine, Z, Gridded(Linear()))
    itp[pts_coarse[1],pts_coarse[2]]
end


end
