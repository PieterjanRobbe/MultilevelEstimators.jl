## report.jl : loads a history file and make a report with the main results

function report(h::History;folder="report $(h[:name])"::AbstractString)

    # problem file name
    fname = isempty(h[:name]) ? "report.tex" : "report_$(h[:name]).tex"

    # make the required directories
    !isdir(folder) && mkdir(folder)
    !isdir(joinpath(folder,"data")) && mkdir(joinpath(folder,"data"))
    !isdir(joinpath(folder,"figures")) && mkdir(joinpath(folder,"figures"))

    # make all data files
    write_E(h,folder)
    write_dE(h,folder)
    write_V(h,folder)
    write_dV(h,folder)
    write_W(h,folder)
    write_samples(h,folder)
    write_index_set(h,folder)
    write_adaptive_index_set(h,folder)
    write_runtime(h,folder)
    write_runtime_vs_rmse(h,folder)
    write_cost(h,folder)

    # make the tikz files
    write_figure_E(h,folder,fname)
    write_figure_V(h,folder,fname)
    write_figure_W(h,folder,fname)
    write_figure_samples(h,folder,fname)
    write_figure_index_set(h,folder,fname)
    write_figure_adaptive_index_set(h,folder,fname)
    write_figure_runtime(h,folder,fname)
    write_figure_runtime_vs_rmse(h,folder,fname)
    write_figure_cost(h,folder,fname)

    # make the report file
    write_main_file(h,folder,fname)
end

for (name,start) in zip(["E","dE","V","dV","W"],[1,2,1,2,1])
    ex = :(
           function $(Symbol("write_",name))(h::History, folder::AbstractString)
               d = h[:ndims]
               for i in 2:2^d
                   proto_idx = ind2sub(tuple([2 for i = 1:d]...),i).-1
                   cntr = 0
                   x = Float64[]
                   y = Float64[]
                   while cntr.*proto_idx ∈ h[:index_set]
                       idx = findfirst([cntr.*proto_idx].==h[:index_set])
                       push!(x,cntr)
                       push!(y,log2(abs(h[Symbol(Symbol($(name)))][idx])))
                       cntr += 1
                   end
                   open(joinpath(folder,"data",string($(name),join(string.(proto_idx)),".txt")), "w") do f
                       for i in $(start):length(x)
                      write(f, @sprintf("%i %12.5f\n",x[i],y[i]))
                  end
              end
          end
      end
      )
    eval(ex)
end

function write_samples(h::History, folder::AbstractString)
    d = h[:ndims]
    if d == 1
        open(joinpath(folder,"data","samples.txt"), "w") do f
            m = maximum(h[:index_set])[1] + 1
            n = h.iter
            S = zeros(m,n)
            for j = 1:n
                s = h[j][:nsamples]
                S[1:length(s),j] = s
            end
            for i in 1:m
                str = @sprintf("%i",i-1)
                for j in 1:n
                    str = string(str,@sprintf(" %i",S[i,j]))
                end
                write(f,string(str,"\n"))
            end
        end
    end
end

function write_index_set(h::History, folder::AbstractString)
    d = h[:ndims]
    tols = [h[i][:tol] for i in 1:h.iter]
    for i in 1:length(tols)
        open(joinpath(folder,"data",string("index_set_",i,".txt")), "w") do f
            index_set = h[i][:index_set]
            for idx in index_set
                for j in 1:d
                    write(f, @sprintf("%i ",idx[j]))
                end
                write(f, "\n")
            end
        end
    end
end

function write_adaptive_index_set(h::History, folder::AbstractString)
    d = h[:ndims]
	dicts = h[:adaptive_index_set]
    for i in 1:length(dicts)
		for (val,name,mode) in zip(0:3,["old","active","active","maximum"],["w" "w" "a" "w"])
			open(joinpath(folder,"data",string("level_",i,"_",name,".txt")), mode) do f
				index_set = collect(keys(dicts[i]))[collect(values(dicts[i])).==val] 
            	for idx in index_set
                	for j in 1:d
                    	write(f, @sprintf("%i ",idx[j]))
                	end
                	write(f, "\n")
            	end
        	end # open
		end
    end
end

function write_runtime_vs_rmse(h::History, folder::AbstractString)
    open(joinpath(folder,"data","runtime_vs_rmse.txt"), "w") do f
        x = [h[i][:rmse] for i in 1:h.iter]
        y = cumsum([h[i][:runtime] for i in 1:h.iter])
        for i in 1:length(x)
            write(f, @sprintf("%12.5f %12.5e\n",x[i],y[i]))
        end
    end
end

function write_runtime(h::History, folder::AbstractString)
    open(joinpath(folder,"data","runtime.txt"), "w") do f
        x = [h[i][:tol] for i in 1:h.iter]
        y = cumsum([h[i][:runtime] for i in 1:h.iter])
        for i in 1:length(x)
            write(f, @sprintf("%12.5f %12.5e\n",x[i],y[i]))
        end
    end
end

function write_cost(h::History, folder::AbstractString)
    open(joinpath(folder,"data","cost.txt"), "w") do f
        x = [h[i][:tol] for i in 1:h.iter]
        y = [h[i][:cost] for i in 1:h.iter]
        for i in 1:length(x)
            write(f, @sprintf("%12.5f %12.5e\n",x[i],y[i]))
        end
    end
end

for (name,color,ylabel) in zip(["E","V"],["red","blue"],["|E[\\;\\cdot\\;]|","V[\\;\\cdot\\;]"])
    ex = :(
           function $(Symbol("write_figure_",name))(h::History, folder::AbstractString, fname::AbstractString)
               open(joinpath(folder,"figures",string($(name),".tex")), "w") do f
                   L = maximum(maximum.(h[:index_set]))
                   d = length(h[:index_set][1])
                   str = tikz_header("level \$\\ell\$",string("\$\\log_2(",$(ylabel),")\$"),
                                     "xmin=0,\n xmax=$(L),\n xtick={0,1,...,$(L)},\n",fname)
                   if d > 1
                       for i in 2:2^d
                           proto_idx = ind2sub(tuple([2 for i = 1:d]...),i).-1
                           marker_idx = mod(i-2,length(PGFPLOTSMARKERS))+1
                           linestyle_idx = mod(sum(proto_idx)-1,length(LINESTYLES))+1
                           name = string($(name),join(string.(proto_idx)))
                           str = string(str,tikz_add_plot($(color),"solid",PGFPLOTSMARKERS[marker_idx],name,"",""))
                           legend_name = string("\$\\Delta Q_{",replace(string(proto_idx),"1","\\ell"),"}\$")
                           str = string(str,tikz_add_plot($(color),LINESTYLES[linestyle_idx],
                                                          PGFPLOTSMARKERS[marker_idx],
                                                          string("d",name),legend_name,"")
                           )
                       end
                   else
                       str = string(str,tikz_add_plot($(color),"solid",PGFPLOTSMARKERS[1],
                                                      string($(name),"1"),"\$Q_\\ell\$",""))
                       str = string(str,tikz_add_plot($(color),LINESTYLES[1],PGFPLOTSMARKERS[1],
                                                      string("d",$(name),"1"),"\$Q_\\ell-Q_\{\\ell-1\}\$",""))
                   end
                   str = string(str,tikz_footer())
                   write(f, str)
               end
           end
           )
    eval(ex)
end

function write_figure_W(h::History, folder::AbstractString, fname::AbstractString)
    open(joinpath(folder,"figures","W.tex"), "w") do f
        L = maximum(maximum.(h[:index_set]))
        d = h[:ndims]
        str = tikz_header("level \$\\ell\$",string("\$\\log_2(C_\\ell)\$"),
                          "xmin=0,\n xmax=$(L),\n xtick={0,1,...,$(L)},\n",fname)
        for i in 2:2^d
            proto_idx = ind2sub(tuple([2 for i = 1:d]...),i).-1
            marker_idx = mod(i-2,length(PGFPLOTSMARKERS))+1
            linestyle_idx = mod(sum(proto_idx)-1,length(LINESTYLES))+1
            name = string("W",join(string.(proto_idx)))
            legend_name = string("\$W_{",replace(string(proto_idx),"1","\\ell"),"}\$")
            str = string(str,tikz_add_plot("green","solid",PGFPLOTSMARKERS[marker_idx],name,legend_name,""))
        end
        str = string(str,tikz_footer())
        write(f, str)
    end
end

function write_figure_samples(h::History, folder::AbstractString, fname::AbstractString)
    d = h[:ndims]
    if d == 1
        open(joinpath(folder,"figures","samples.tex"), "w") do f
            L = maximum(maximum.(h[:index_set]))
            str = tikz_header("level \$\\ell\$","number of samples \$N_\\ell\$",
                              "ymode=log,\n xmin=0,\n xmax=$(L),\n xtick={0,1,...,$(L)},\n",fname)
            cm = ColorMap("jet")
            v = linspace(0,1,h.iter)
            for i in 1:h.iter
                color = cm(v[i])
                color_name = "{rgb,1:red,$(color[1]); green,$(color[2]); blue,$(color[3])}"
                str = string(str,tikz_add_plot(color_name,"solid","*","samples","","x index = {0}, y index = {$(i)}"))
            end
            str = string(str,tikz_footer())
            write(f, str)
        end
    end
end

function write_figure_index_set(h::History, folder::AbstractString, fname::AbstractString)
    d = h[:ndims]
    if d == 2
        write_figure_index_set_2(h, folder, fname)
    elseif d == 3
        write_figure_index_set_3(h, folder, fname)
    end
end

for d in ["2" "3"]
    ex = :(function $(Symbol("write_figure_index_set_",d))(h::History, folder::AbstractString, fname::AbstractString)
               for i = 1:h.iter
                   open(joinpath(folder,"figures","index_set_$(i).tex"), "w") do f
                       max_level = maximum(maximum.(h[:index_set]))
                       str = $(Symbol("tikz_index_set_",d,"d"))(fname,i,max_level,false,false,false,false,1,"")
                       write(f, str)
                   end
               end
               open(joinpath(folder,"figures","index_set.tex"), "w") do f
                   tols = [h[i][:tol] for i in 1:h.iter]
                   str = tikz_index_set_table(fname,tols)
                   write(f, str)
               end
           end
          )
    eval(ex)
end

function write_figure_adaptive_index_set(h::History, folder::AbstractString, fname::AbstractString)
    d = h[:ndims]
	if isa(h[:method],AD)
    	if d == 2
        	write_figure_adaptive_index_set_2(h, folder, fname)
    	elseif d == 3
        	write_figure_adaptive_index_set_3(h, folder, fname)
    	end
	end
end

for d in ["2" "3"]
	ex = :(function $(Symbol("write_figure_adaptive_index_set_",d))(h::History, folder::AbstractString, fname::AbstractString)
			   dicts = h[:adaptive_index_set]
			   for i = 1:length(dicts)
				   open(joinpath(folder,"figures","adaptive_index_set_$(i).tex"), "w") do f
					   max_level = maximum(maximum.(keys(dicts[end])))
					   v = collect(values(dicts[i]))
					   has_active = any(v.==1) || any(v.==2)
					   has_old = any(v.==0)
					   has_maximum = any(v.==3)
					   str = $(Symbol("tikz_index_set_",d,"d"))(fname,i,max_level,true,has_active,has_old,has_maximum,1,"")
					   write(f, str)
				   end
			   end
			   open(joinpath(folder,"figures","adaptive_index_set.tex"), "w") do f
				   str = tikz_adaptive_index_set_table(fname,0:length(dicts)-1)
				   write(f, str)
			   end
			   for mode in ["active","old","maximum"]
				   open(joinpath(folder,"data","index_set_legend_$(mode).txt"), "w") do f
					   dim = h[:ndims]
					   for j in 1:dim
						   write(f, @sprintf("%i ",0))
					   end
				   end
				   open(joinpath(folder,"figures","index_set_legend_$(mode).tex"), "w") do f
					   max_level = maximum(maximum.(keys(dicts[end])))
					   scaling = (max_level+1.2)/1.2
					   write(f, $(Symbol("tikz_index_set_",d,"d"))(fname,"legend",0,false,false,false,false,scaling,mode))
				   end
			   end
		   end
		  )
	eval(ex)
end

for (name,color,ylabel) in zip(["runtime","cost","runtime_vs_rmse"],["blue","red","blue"],["run time","standard cost","run time"])
    ex = :(
           function $(Symbol("write_figure_",name))(h::History, folder::AbstractString, fname::AbstractString)
               open(joinpath(folder,"figures",string($(name),".tex")), "w") do f
                   str = tikz_header("accuracy",$(ylabel),"xmode = log,\n ymode = log",fname)
                   str = string(str,tikz_add_plot($(color),"solid","*",$(name),"",""))
                   str = string(str,tikz_footer())
                   write(f, str)
               end
           end
           )
    eval(ex)
end

function write_main_file(h::History,folder::AbstractString,fname::AbstractString)
    open(joinpath(folder,fname), "w") do f
		write(f, file_contents(h[:name],length(h[:index_set][1]),isa(h[:method],AD)))
    end
end
