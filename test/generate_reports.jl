## automatically generate all reports from the test files

using JLD, MultilevelEstimators, PyPlot, Suppressor

for folder in readdir("data")
    folderpath = joinpath("data",folder)
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
