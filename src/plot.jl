## plot.jl : loads a history file and make some plots

function plot_E(h::History)
    L = h[:max_level]
    x1 = 0:L
    y1 = log2.(abs.(h[:E]))
    x2 = 1:L
    y2 = log2.(abs.(h[:dE]))[2:end]
    figure()
    plot(x1,y1,linestyle="-",marker="o",color="blue",label="E")
    plot(x2,y2,linestyle="--",marker="o",color="blue",label="ΔE")
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
    L = h[:max_level]
    x1 = 0:L
    y1 = log2.(h[:V])
    x2 = 1:L
    y2 = log2.(h[:dV])[2:end]
    figure()
    plot(x1,y1,linestyle="-",marker="o",color="red",label="V")
    plot(x2,y2,linestyle="--",marker="o",color="red",label="ΔV")
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
    L = h[:max_level]
    x = 0:L
    y = log2.(h[:W])
    figure()
    plot(x,y,linestyle="-",marker="o",color="green")
    ax = gca()
    m = matplotlib[:ticker][:MultipleLocator](1)
    ax[:xaxis][:set_major_locator](m)
    xlabel("level")
    ylabel("log₂(⋅)")
    title("increase in W")
    show()
end

function plot_samples(h::History)
    L = h[:max_level]
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
