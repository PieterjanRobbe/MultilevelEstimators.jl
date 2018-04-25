## automatically generate all reports from the test files

using JLD, MultilevelEstimators, PyPlot, Suppressor

data_dir = string(Pkg.dir("MultilevelEstimators"),"applications/SPDE/data")
for folder in readdir(data_dir)
    folderpath = joinpath(data_dir,folder)
    filepath = joinpath(folderpath,"history.jld")
    if isfile(filepath)
        print(string("generating report for ",folder,"..."))
        h = load(filepath,"history"); 
        report(h,folder=folderpath)
        println("done")
        cd(folderpath)
        texfile = string("report_",folder,".tex")
        print(string("running pdflatex..."))
        @suppress run(`pdflatex $(texfile)`)
        println("done")
        cd("../..")
    end
end
