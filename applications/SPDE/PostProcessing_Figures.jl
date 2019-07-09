module PostProcessing_Figures

using   Base.Test, JLD,PyPlot, Distributions


data_dir = string(joinpath(Pkg.dir("MultilevelEstimators"),"applications","SPDE","data"))

#for folder in readdir(data_dir)
folder="Beam MLMC (multiple)"
    folderpath = joinpath(data_dir,folder)
    filepath = joinpath(folderpath,"history.jld")

    println(filepath)
    if isfile(filepath)
#        print(string("generating report for ",folder,"..."))
        h = load(filepath,"history");
        nb_qoi=size(h[:samples])[1]
        nb_of_levels=size(h[:index_set])[1]
        println(nb_qoi)
        println(nb_of_levels)
        println(h[:E][1])
        println(h[:E][2])
        println(h[:E][3])
        println(h[:E][4])

        AverageValueDisplacementY=zeros(nb_qoi,1)


        for idx=1:nb_qoi

        for id=0:nb_of_levels-1

        AverageValueDisplacementY[idx,1]=AverageValueDisplacementY[idx,1]+(sum(h[:samples][idx][(id,)]))/size(h[:samples][idx][(id,)])[1]

        end
        end

        pdf_vec=zeros(1000,nb_qoi)

        for outloop=1:1000

        for idx=1:nb_qoi

        for id=0:nb_of_levels-1

        pdf_vec[outloop,idx]=rand((h[:samples][idx][(id,)]))+pdf_vec[outloop,idx]
        #println(rand((h[:samples][idx][(id,)])))

        end
        end
        end
      for id=2:nb_qoi-1

      pdfSample=pdf_vec[1:end,id]

    #   println(pdfSample)
       expo=ceil(abs(log(minimum((pdfSample)))/log(10)))+4;
       println(expo)
       fittedDist=Distributions.fit(Normal,pdfSample)
    #   println(fittedDist)

       ValMax=abs(maximum(pdfSample));
       ValMin=abs(minimum(pdfSample));

      if(ValMax>ValMin)
      Val=ValMax;
      else
     Val=ValMin;
       end
      divider=10
   #   println(minimum(pdfSample)-Val/divider)
   #   println(5^-expo)
   #   println(maximum(pdfSample)+Val/divider)
      startpoint=minimum(pdfSample)-Val/divider
      endpoint=maximum(pdfSample)+Val/divider
      x_values=startpoint:(5.0^-expo):endpoint;

      y_values=pdf.(fittedDist,x_values)
      #figure()
      #plot(x_values,y_values)
      end

#fill_between(vec(collect(x)),vec(collect(y)),vec(-collect(y)),facecolor="blue", alpha=0.1)

        figure()
        plot(vec(collect(1:1:nb_qoi)),-AverageValueDisplacementY)

        #println((sum(h[:samples][81][(0,)]))/size(h[:samples][81][(0,)])[1])
        #println((sum(h[:samples][81][(1,)]))/size(h[:samples][81][(1,)])[1])
        #println((sum(h[:samples][81][(2,)]))/size(h[:samples][81][(2,)])[1])
        #println((sum(h[:samples][81][(3,)]))/size(h[:samples][81][(3,)])[1])


        #println(h[:dE][2])
        #println(h[:dE][3])
        #println(h[:dE][4])


        #println(size(h[:samples][80][(0,)]))
        #println(size(h[:index_set]))
        #println(h[:W])
        #println(size(h[:samples]))


    end
#end





end
