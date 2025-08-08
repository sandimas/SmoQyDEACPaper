using SmoQyDEAC
using FileIO

if Threads.nthreads() == 1
    println("Rerun this with --threads=auto")
end


function run_deac(folder)
    n_bin = 100
    n_run = 20
    dict = load(joinpath(folder,"true.jld2"))
    G = dict["Gτ"]
    β = dict["β"]
    τs = collect(dict["τs"])
    ωs = collect(LinRange(-10.0,10.0,400))
    ferm = DEAC_Binned(G,β,τs,ωs,"time_fermionic",n_bin,n_run,joinpath(folder,"deac.jld2"),joinpath("distributions",folder,"chk.jld2"))
end



folders  = ["01_Single_Gauss","02_Single_Lorentz","03_Double_Gauss","04_Double_Lorentz","05_Double_mixed",
            "06_Double_Plateau","07_Triple_1","08_Triple_2","09_Chaos","10_Batman"];

for folder in folders
    run_deac(folder)
end
