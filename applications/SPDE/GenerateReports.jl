#### COMMENT : Shut down Julia , run the piece of code in comment , then execute file
module GenerateReports
using Pkg
using Reporter
using FileIO

data_dir = string(joinpath(Pkg.dir("MultilevelEstimators"),"applications","SPDE","data"))
filepath = joinpath(data_dir,"Beam MC_Hom_Lin_Ref")
filepath_2 = string(filepath,"/Beam MC_Hom_Lin_Ref.jld2")
history=load(filepath_2,"history")
report(history,filepath,include_preamble=true)

data_dir = string(joinpath(Pkg.dir("MultilevelEstimators"),"applications","SPDE","data"))
for folder in readdir(data_dir)
    folderpath = joinpath(data_dir,folder)
     endname=string(basename(folderpath))q
    filepath = string(folderpath,"/")
    filepath=string(filepath,endname)
    filepath_file = string(filepath,".jld2")
    println(filepath_file)

    println(filepath_file)
    println(isfile(filepath_file))
    if isfile(filepath_file)
        history = load(filepath_file,"history");
        report(history,filepath)

#    #    println("done")
##        cd(folderpath)
# #       texfile = string("report_",folder,".tex")
#  #      print(string("running pdflatex..."))
#  #      @suppress run(`pdflatex $(texfile)`)
#  #      println("done")
#  #      cd("../..")
    end
end
end
