using ACFlow
using FileIO
using Statistics

function run_mem(folder)
    dict = load(joinpath(folder,"true.jld2"))
    
    
    G = mean(dict["G"],dims=1)[1,:]
    err = std(dict["G"],dims=1)[1,:]
    β = dict["β"]
    τs = dict["τs"]
    nτ = size(τs,1)
    ωmax = dict["ωs"][end]
    
    B_dict = Dict{String,Any}(
        "solver" => "MaxEnt",
        "ktype"  => "bsymm",
        "mtype"  => "flat",
        "mesh"   => "linear",
        "grid"   => "btime",
        "nmesh"  => 400,
        "ngrid"  => nτ,
        "wmax"   => ωmax,
        "wmin"   => -ωmax,
        "beta"   => β,  
        "fwrite" => false
    )
    S_dict = Dict{String,Any}( )
    ACFlow.setup_param(B_dict,S_dict)
    
    ω_MEM, A_MEM, _ = ACFlow.solve(collect(τs),G,err .+ 0.00000001)

    mem = Dict{String,Any}( "A"=>A_MEM, "ωs"=>ω_MEM)

    save(joinpath(folder,"mem.jld2"),mem)
end

folders = ["01_Gauss_1.0","02_Gauss_0.1","03_Gauss_0.01",
            "04_Lorentz_1.0","05_Double_Gauss_1","06_Double_Gauss_2",
            "07_Gauss_Singular","08_Log_Singular","09_Chaos","10_Chaos_Singular"]

for folder in folders
    run_mem(folder)
end