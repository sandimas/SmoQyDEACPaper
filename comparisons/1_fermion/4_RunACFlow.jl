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
    
    B_dict = Dict{String,Any}(
        "solver" => "MaxEnt",
        "ktype"  => "fermi",
        "mtype"  => "flat",
        "mesh"   => "linear",
        "grid"   => "ftime",
        "nmesh"  => 400,
        "ngrid"  => nτ,
        "wmax"   => 10.0,
        "wmin"   => -10.0,
        "beta"   => β,  
        "fwrite" => false
    )
    S_dict = Dict{String,Any}( )
    ACFlow.setup_param(B_dict,S_dict)
    
    ω_MEM, A_MEM, _ = ACFlow.solve(collect(τs),G,err .+ 0.00000001)

    mem = Dict{String,Any}( "A"=>A_MEM, "ωs"=>ω_MEM)

    save(joinpath(folder,"mem.jld2"),mem)
end


folders  = ["01_Single_Gauss","02_Single_Lorentz","03_Double_Gauss","04_Double_Lorentz","05_Double_mixed",
            "06_Double_Plateau","07_Triple_1","08_Triple_2","09_Chaos","10_Batman"]

for folder in folders
    run_mem(folder)
end