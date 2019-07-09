 module FieldTransformation

 #using Distributions,PyPlot,Interpolations
 using Distributions,Interpolations, GaussianRandomFields



function Transform_NL(Points::Array{Float64,2})
Cov=0.0327
beta=200E3*Cov^2/(1-Cov^2);
alpha=200E3/beta;
GamDist=Gamma(alpha,beta)
#@show Cov
Evaluation_points_Gam=(1.6:0.001:2.34)*10.0E4;
CDF_Gamma=cdf.(GamDist,Evaluation_points_Gam);
NormDist=Normal(0.0,1.0)
Evaluation_points_Norm=-4.0:5.0^-4:4.0;
CDF_Norm=cdf.(NormDist,Evaluation_points_Norm)
cdf_mat=interp1(Evaluation_points_Norm,CDF_Norm,Points)
GammaField=interp1(CDF_Gamma,Evaluation_points_Gam,cdf_mat)

#@show CDF_Gamma
#figure()
#plot(Evaluation_points_Gam,CDF_Gamma)

#figure()
#plot(Evaluation_points_Norm,CDF_Norm)

return GammaField
end

function Transform_L(Points::Array{Float64,2})
Cov=0.35
beta=30E3*Cov^2/(1-Cov^2);
alpha=30E3/beta;
GamDist=Gamma(alpha,beta)
#@show Cov
Evaluation_points_Gam=(0:0.01:10)*10.0E4;
CDF_Gamma=cdf.(GamDist,Evaluation_points_Gam);
NormDist=Normal(0.0,1.0)
Evaluation_points_Norm=-4.0:5.0^-4:4.0;
CDF_Norm=cdf.(NormDist,Evaluation_points_Norm)
cdf_mat=interp1(Evaluation_points_Norm,CDF_Norm,Points)
GammaField=interp1(CDF_Gamma,Evaluation_points_Gam,cdf_mat)

#@show CDF_Gamma
#figure()
#plot(Evaluation_points_Gam,CDF_Gamma)

#figure()
#plot(Evaluation_points_Norm,CDF_Norm)

return GammaField
end



function interp1(xpt, ypt, x)
        y = zeros(size(x,1),size(x,2))
        idx = trues(size(x,1),size(x,2))


        intf = extrapolate(interpolate((xpt,), ypt, Gridded(Linear())),Flat())

         y[idx] = intf(x[idx])



    return y
end

function interp1_equi(xpt, ypt, div; method="linear", extrapvalue=nothing)
    xpt=vec(xpt)
    inc=(xpt[2]-xpt[1])/div
    x=xpt[1]:inc:xpt[size(xpt)[1]]
    ypt=vec(ypt)

    if extrapvalue == nothing
        y = zeros(x)
        idx = trues(x)
    else
        y = extrapvalue*ones(x)
        idx = (x .>= xpt[1]) .& (x .<= xpt[end])
    end

    if method == "linear"
        intf = interpolate((xpt,), ypt, Gridded(Linear()))
        y[idx] = intf[x[idx]]

    elseif method == "cubic"
        itp = interpolate(ypt, BSpline(Cubic(Natural())), OnGrid())
        intf = Interpolations.scale(itp, xpt)
        y[idx] = [intf[xi] for xi in x[idx]]
    end

    return y
end

function Test()

    Lx = 2500.0;
    Ly = 250.0;
    he = 250/32;

    nelx=Lx/he
    nely=Ly/he
    exp_field = GaussianRandomFields.Exponential(0.3,σ=1.0,p=1.0)


    cov = GaussianRandomFields.CovarianceFunction(2,exp_field)
    #srand(1234)
    grf = GaussianRandomFields.GaussianRandomField(cov,KarhunenLoeve(101),0:he/Lx:0.99999,0:he/Ly:0.99999,quad=GaussLegendre())

    #contourf(grf)
    #plot(grf)
    #@show GaussianRandomFields.sample(grf)

    E = 200E3;
    nu = 0.25;
    fy = 240.0;
    t = 1.0;
#    srand(1234)
    #@show ra
    #GammaDist=Gamma(16,16)
    #RandomGammaNumber=rand(GammaDist,1)
    #@show RandomGammaNumber
    for id=1:10
        ra=randn(101)

   f1= GaussianRandomFields.sample(grf;xi=ra);
#    figure()
#    surf(f1)

    E=f1
    GammaField=FieldTransformation.Transform_NL(f1);
#    figure()
#    surf(GammaField)

    Zc = Array(view(GammaField, 2:2:size(GammaField, 1), 2:2:size(GammaField, 2)))

#    figure()
#    surf(Zc)
end
#Check that GammaField does not contain NAN?
#    Dispx,Dispy,u = mxcall(:Solver_NL_JULIA_MATLAB,3,Lx,Ly,he,GammaField,nu,fy,t)


end

function DropValues(InArray)
id=1
id_dub=1
sizeAr=Int64(round((size(InArray)[1])/2))+1
#println(sizeAr)
Result=zeros(sizeAr,1)

while(id_dub<=sizeAr)
Result[id_dub]=InArray[id]
id_dub=id_dub+1
id=id+2

end
return vec(Result)
end

end
