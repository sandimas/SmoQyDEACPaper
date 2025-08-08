using FileIO
using CSV
using DataFrames

folders  = ["01_Single_Gauss","02_Single_Lorentz","03_Double_Gauss","04_Double_Lorentz","05_Double_mixed",
            "06_Double_Plateau","07_Triple_1","08_Triple_2","09_Chaos","10_Batman"];

function maxent(folder)
    dict = load(joinpath(folder,"true.jld2"))
    G = dict["G"]
    N_bin = size(G,1)
    τs = collect(dict["τs"])

    fname = joinpath(folder,"maxent_input.txt")
    open(fname,"w") do file
        for (τi, τ) in enumerate(τs)
            print(file,τ)
            for bin in 1:N_bin
                print(file," ",G[bin,τi])
            end
            print(file,"\n")
        end
    end

    mycommand = `python3 ../9_scripts/Maxent_b.py $fname $folder/`
    
    run(mycommand)

    df = CSV.read(joinpath(folder,"maxent.dat"),DataFrame;header=false)

    ωs = df[:,1]
    A = df[:,2] 
    outdict = Dict{String,Any}("ωs"=>ωs,"A"=>A)
    save(joinpath(folder,"mem_cov.jld2"),outdict)

end

for folder in folders

    maxent(folder)
end
