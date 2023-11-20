using SynthAC
using FileIO

folders = ["01_Gauss_1.0","02_Gauss_0.1","03_Gauss_0.01",
            "04_Lorentz_1.0","05_Double_Gauss_1","06_Double_Gauss_2",
            "07_Gauss_Singular","08_Log_Singular","09_Chaos","10_Chaos_Singular"]


for i in 1:size(folders,1)
    folders[i] = joinpath("distributions",folders[i])
    try
        mkdir(folders[i])
    catch
    end
end

#
folder = folders[1]
dists = []

AppendDistribution!(dists,Normal(1.0,0.2,A=0.5))

dict = GenerateCorrelationFunctions(dists,10.0,0.05,false;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000,Maxω=3.0)

#
folder = folders[2]
dists = []

AppendDistribution!(dists,Normal(0.1,0.05,A=0.75))

dict = GenerateCorrelationFunctions(dists,10.0,0.05,false;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000,Maxω=1.0)

#
folder = folders[3]
dists = []

AppendDistribution!(dists,Normal(0.01,0.05,A=0.75))

dict = GenerateCorrelationFunctions(dists,10.0,0.05,false;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000,Maxω=1.0)

#
folder = folders[4]
dists = []

AppendDistribution!(dists,Cauchy(1.0,0.5,A=0.75))

dict = GenerateCorrelationFunctions(dists,10.0,0.05,false;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000,Maxω=10.0)

#
folder = folders[5]
dists = []

AppendDistribution!(dists,Normal(1.0,0.3,A=0.5))
AppendDistribution!(dists,Normal(2.5,0.4,A=0.25))

dict = GenerateCorrelationFunctions(dists,10.0,0.05,false;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000,Maxω=5.0)

#
folder = folders[6]
dists = []

AppendDistribution!(dists,Normal(1.40,0.3,A=0.5))
AppendDistribution!(dists,Normal(4.,1.,A=0.25))

dict = GenerateCorrelationFunctions(dists,10.0,0.05,false;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000,Maxω=8.0)

#
folder = folders[7]
dists = []

function gauss_sing(x)
    if x < 0.5
        return 0.0
    else
       return (0.01/(abs(x)-0.49999) +  Normal(2.,0.5,A=1.0)(x))/abs(x)
    end
end

AppendDistribution!(dists,gauss_sing)

dict = GenerateCorrelationFunctions(dists,10.0,0.05,false;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000,Maxω=4.0)

#
folder = folders[8]
dists = []

function gauss_sing2(x)
    if abs(x) < 0.5
        return 0.0
    elseif abs(x) > 3.0
        return 0.0
    else
       return (0.5/(abs(x)-0.499) )*(3.0-abs(x))/(abs(x))
    end
end

AppendDistribution!(dists,gauss_sing2)


dict = GenerateCorrelationFunctions(dists,10.0,0.05,false;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000,Maxω=4.0)

#
folder = folders[9]
dists = []


AppendDistribution!(dists,Normal(1.40,0.3,A=0.5))
AppendDistribution!(dists,Normal(2.,1.,A=0.25))
AppendDistribution!(dists,x->Normal(3.,1.,A=0.25)(x) * (sin(5.5*x)^2 + abs(cos(3.0*x))))


dict = GenerateCorrelationFunctions(dists,10.0,0.05,false;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000,Maxω=10.0)

#
folder = folders[10]
dists = []

function wild(x)
    if x < 0.5
        return 0
    else
        return (1 /(x-0.4999)) * (Normal(1.40,0.3,A=0.5)(x)+Normal(2.,1.,A=0.25)(x)+Normal(3.,1.,A=1.0)(x) * (sin(5.5*x)^2 + abs(cos(3.0*x))))
    end
end

AppendDistribution!(dists,x->wild(x))

dict = GenerateCorrelationFunctions(dists,10.0,0.05,false;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000,Maxω=10.0)