using FileIO
using Printf
using DataFrames
using CSV


# To get the SAC code I used download and build from here
# https://github.com/JefferyWangSH/sac.git


folders  = ["01_Single_Gauss","02_Single_Lorentz","03_Double_Gauss","04_Double_Lorentz","05_Double_mixed",
            "06_Double_Plateau","07_Triple_1","08_Triple_2","09_Chaos","10_Batman"];



function sac(folder)
    
    cp("sac_files/config.toml",joinpath(folder,"config.toml"))
    dict = load(joinpath(folder,"true.jld2"))
    G = dict["G"]
    N_bin = size(G,1)
    τs = collect(dict["τs"])

    open(joinpath(folder,"tgrids.in"),"w") do file
        println(file,@sprintf("%20u%20f",size(τs,1),τs[end]))
        for (τi, τ) in enumerate(τs)
            println(file,@sprintf("%20u%20f",τi-1,τ))
        end
    end
    open(joinpath(folder,"corr.in"),"w") do file
        println(file,@sprintf("%20u%20u%20f",N_bin,size(τs,1),τs[end]))
        for bin in 1:N_bin
            for (τi, τ) in enumerate(τs)
                println(file,@sprintf("%20u%20u%20f",bin-1,τi-1,G[bin,τi]))
            end
        end
    end

    config = joinpath(folder,"config.toml")
    tgrids = joinpath(folder,"tgrids.in")
    corr = joinpath(folder,"corr.in")
    log =  joinpath(folder,"sac.log")
    spec = joinpath(folder,"sac.out")
    report = joinpath(folder,"sac.rpt")

    command = `sac --config=$config --tgrids=$tgrids --corr=$corr --log=$log --spec=$spec --report=$report`
    
    run(command)

 
    df = CSV.read(joinpath(folder,"sac.out"),DataFrame;delim=" ",ignorerepeated=true,header=false)

    ωs = df[:,2]
    A = df[:,3] 

    outdict = Dict{String,Any}("ωs"=>ωs,"A"=>A)
    save(joinpath(folder,"sac.jld2"),outdict)

end

for folder in folders
    sac(joinpath("distributions",folder))
end

