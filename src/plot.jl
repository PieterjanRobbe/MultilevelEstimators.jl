## plot.jl : loads a history file and make some plots

const MARKERS = ["o","s","v","^","<",">","*"]
const COLORS = ["blue","green","red"]

function plot_E(h::History)
	d = length(h[:index_set][1])
	figure()
	for i in 2:2^d
		proto_idx = ind2sub(tuple([2 for i = 1:d]...),i).-1
		cntr = 0
		x1 = Float64[]; x2 = Float64[]
		y1 = Float64[]; y2 = Float64[]
		while cntr.*proto_idx ∈ h[:index_set]
			idx = findfirst([cntr.*proto_idx].==h[:index_set])
			push!(x1,cntr); push!(x2,cntr)
			push!(y1,log2(abs(h[:E][idx]))),push!(y2,log2(abs(h[:dE][idx])))
			cntr += 1
		end
		plot(x1,y1,linestyle="-",marker=MARKERS[mod(i-2,length(MARKERS))+1],color=COLORS[mod(sum(proto_idx)-1,length(COLORS))+1],label="E, $(replace(string(proto_idx),"1","\u2605"))")
		plot(x2[2:end],y2[2:end],linestyle="--",marker=MARKERS[mod(i-2,length(MARKERS))+1],color=COLORS[mod(sum(proto_idx)-1,length(COLORS))+1],label="ΔE, $(replace(string(proto_idx),"1","\u2605"))")
	end
	ax = gca()
	m = matplotlib[:ticker][:MultipleLocator](1)
	ax[:xaxis][:set_major_locator](m)
	xlabel("level")
	ylabel("log₂(|⋅|)")
	legend(loc="lower left",fancybox="true")
	title("decay of E")
	show()
end

function plot_V(h::History)
	d = length(h[:index_set][1])
	figure()
	for i in 2:2^d
		proto_idx = ind2sub(tuple([2 for i = 1:d]...),i).-1
		cntr = 0
		x1 = Float64[]; x2 = Float64[]
		y1 = Float64[]; y2 = Float64[]
		while cntr.*proto_idx ∈ h[:index_set]
			idx = findfirst([cntr.*proto_idx].==h[:index_set])
			push!(x1,cntr); push!(x2,cntr)
			push!(y1,log2(h[:V][idx])),push!(y2,log2(h[:dV][idx]))
			cntr += 1
		end
		plot(x1,y1,linestyle="-",marker=MARKERS[mod(i-2,length(MARKERS))+1],color=COLORS[mod(length(COLORS)-sum(proto_idx),length(COLORS))+1],label="V, $(replace(string(proto_idx),"1","\u2605"))")
		plot(x2[2:end],y2[2:end],linestyle="--",marker=MARKERS[mod(i-2,length(MARKERS))+1],color=COLORS[mod(length(COLORS)-sum(proto_idx),length(COLORS))+1],label="ΔV, $(replace(string(proto_idx),"1","\u2605"))")
	end
	ax = gca()
	m = matplotlib[:ticker][:MultipleLocator](1)
	ax[:xaxis][:set_major_locator](m)
	xlabel("level")
	ylabel("log₂(⋅)")
	legend(loc="lower left",fancybox="true")
	title("decay of V")
	show()
end

function plot_W(h::History)
	d = length(h[:index_set][1])
	figure()
	for i in 2:2^d
		proto_idx = ind2sub(tuple([2 for i = 1:d]...),i).-1
		cntr = 0
		x = Float64[]
		y = Float64[]
		while cntr.*proto_idx ∈ h[:index_set]
			idx = findfirst([cntr.*proto_idx].==h[:index_set])
			push!(x,cntr)
			push!(y,log2(h[:W][idx]))
			cntr += 1
		end
		plot(x,y,linestyle="-",marker=MARKERS[mod(i-2,length(MARKERS))+1],color=COLORS[mod(sum(proto_idx)+1,length(COLORS))+1],label="W, $(replace(string(proto_idx),"1","\u2605"))")
	end
	ax = gca()
	m = matplotlib[:ticker][:MultipleLocator](1)
	ax[:xaxis][:set_major_locator](m)
	xlabel("level")
	ylabel("log₂(⋅)")
	legend(loc="lower left",fancybox="true")
	title("increase in W")
	show()
end

function plot_samples(h::History)
	d = length(h[:index_set][1])
	if d == 1
		fig = figure()
		cm = ColorMap("jet")
		v = linspace(0,1,h.iter)
		for iter = 1:h.iter
			y = h[iter][:nsamples]
			x = 0:length(y)-1
			plot(x,y,linestyle="-",marker="o",color=cm(v[iter]))
		end
		ax = gca()
		m = matplotlib[:ticker][:MultipleLocator](1)
		ax[:xaxis][:set_major_locator](m)
		ax[:set_yscale]("log")
		xlabel("level")
		ylabel("number of samples")
		legend(loc="lower left",fancybox="true")
		title("number of samples")
		show()
	end
end

function plot_time(h::History)
    x = [h[i][:tol] for i in 1:h.iter]
    y = cumsum([h[i][:runtime] for i in 1:h.iter])
    figure()
    plot(x,y,linestyle="-",marker="o",color="blue")
    ax = gca()
    ax[:set_xscale]("log")
    ax[:set_yscale]("log")
    xlabel("tolerance")
    ylabel("run time")
    title("run time")
    show()
end

function plot_cost(h::History)
    x = [h[i][:tol] for i in 1:h.iter]
    y = [h[i][:cost] for i in 1:h.iter]
    figure()
    plot(x,y,linestyle="-",marker="o",color="red")
    ax = gca()
    ax[:set_xscale]("log")
    ax[:set_yscale]("log")
    xlabel("tolerance")
    ylabel("standard cost")
    title("standard cost")
    show()
end
