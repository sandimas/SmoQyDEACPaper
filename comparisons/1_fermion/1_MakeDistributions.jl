try 
    using SynthAC
catch
    using Pkg
    Pkg.add(url="https://github.com/sandimas/SynthAC.jl.git")
    using SynthAC
end

using FileIO

#
folder = "01_Single_Gauss"
mkpath(folder)

dists = []

AppendDistribution!(dists,Normal(1.5,0.5,A=0.75))

dict = GenerateCorrelationFunctions(dists,10.0,0.05,true;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000)

#
folder = "02_Single_Lorentz"
mkpath(folder)
dists = []

AppendDistribution!(dists,Cauchy(-1.5,0.5,A=0.75))

dict = GenerateCorrelationFunctions(dists,10.0,0.05,true;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000)

#
folder = "03_Double_Gauss"
mkpath(folder)
dists = []

AppendDistribution!(dists,Normal(-3.5,0.5,A=0.75))
AppendDistribution!(dists,Normal(1.5,1.0,A=0.5))

dict = GenerateCorrelationFunctions(dists,10.0,0.05,true;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000)

#
folder = "04_Double_Lorentz"
mkpath(folder)
dists = []

AppendDistribution!(dists,Cauchy(-3.5,0.5,A=0.35))
AppendDistribution!(dists,Cauchy(1.5,1.0,A=0.75))

dict = GenerateCorrelationFunctions(dists,10.0,0.05,true;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000)

#
folder = "05_Double_mixed"
mkpath(folder)
dists = []

AppendDistribution!(dists,Normal(-4.5,1.0,A=0.35))
AppendDistribution!(dists,Cauchy(3.5,0.5,A=0.75))

dict = GenerateCorrelationFunctions(dists,10.0,0.05,true;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000)

#
folder = "06_Double_Plateau"
mkpath(folder)
dists = []

AppendDistribution!(dists,x->min(1.0,0.2/(x-4)^2))
AppendDistribution!(dists,x->min(1.0,0.2/(x+4)^2))

dict = GenerateCorrelationFunctions(dists,10.0,0.05,true;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000)

#
folder = "07_Triple_1"
mkpath(folder)
dists = []

AppendDistribution!(dists,Normal(-3.5,0.5,A=0.75))
AppendDistribution!(dists,Normal(1.5,1.0,A=0.5))
AppendDistribution!(dists,Normal(4.5,0.6,A=1))

dict = GenerateCorrelationFunctions(dists,10.0,0.05,true;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000)

#
folder = "08_Triple_2"
mkpath(folder)
dists = []

AppendDistribution!(dists,Cauchy(-4.5,0.5,A=0.45))
AppendDistribution!(dists,Normal(0.5,0.33,A=0.75))
AppendDistribution!(dists,Normal(4.5,0.6,A=1))

dict = GenerateCorrelationFunctions(dists,10.0,0.05,true;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000)

#
folder = "09_Chaos"
mkpath(folder)
mkpath(folder)
dists = []

AppendDistribution!(dists,Cauchy(-4.5,0.5,A=0.45))
AppendDistribution!(dists,Normal(0.5,0.33,A=0.75))
AppendDistribution!(dists,Normal(4.5,0.6,A=1))
AppendDistribution!(dists,x->min(1.0,0.2/(x-2.3)^2))
AppendDistribution!(dists,x->min(1.0,0.2/(x+4)^2))
AppendDistribution!(dists,x->min(1.0,1/(x)^2))              

dict = GenerateCorrelationFunctions(dists,10.0,0.05,true;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000)

#
folder = "10_Batman"
mkpath(folder)
dists = []

AppendDistribution!(dists,x->min(1.0,1.0/(x)^2))
AppendDistribution!(dists,Normal(1.5,0.25,A=0.75))
AppendDistribution!(dists,Normal(-1.5,0.25,A=0.75))      

dict = GenerateCorrelationFunctions(dists,10.0,0.05,true;outfile=folder*"/true.jld2",σ0 =0.01,NBins=1000)