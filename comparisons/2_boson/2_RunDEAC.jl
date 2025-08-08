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
    ferm = DEAC_Binned(G,β,τs,ωs,"time_bosonic_symmetric",n_bin,n_run,joinpath(folder,"deac.jld2"),joinpath("distributions",folder,"chk.jld2"))
end


folders = ["01_Gauss_1.0","02_Gauss_0.1","03_Gauss_0.01",
            "04_Lorentz_1.0","05_Double_Gauss_1","06_Double_Gauss_2",
            "07_Gauss_Singular","08_Log_Singular","09_Chaos","10_Chaos_Singular"]


for folder in folders
    run_deac(folder)
end