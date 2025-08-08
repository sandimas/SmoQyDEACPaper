using FileIO
using SmoQyDEAC
include("../9_scripts/SQMCloader.jl")

directory = "new_1D_w1.00_l0.50_n1.00_L32_b40.00-1/"
if Threads.nthreads() == 1
    println("This should be run \$ julia --threads=auto 2_RunDEAC.jl")
end

kpt = 1:32

ωmin = -10.0
ωmax = 10.0
nω = 200
ωs = collect(LinRange(ωmin,ωmax,nω))

outdir = "data/DEAC/"

mkpath(outdir)

if isfile("data/dqmc_data.jld2")
    input_dict = load("data/dqmc_data.jld2")
    G = input_dict["G"]
    τs = input_dict["τs"]
    β = input_dict["β"]
else

    dict,real,imag,_,_,β = load_from_SmoQyDQMC(
        simulationfolder=directory,
        correlation="greens",
        space="momentum",
        type="time_displaced",bin=true
    )

    G = real[1,1,:,:,1,1,:,1]
    τs=collect(LinRange(0.0,β,size(G,1)))
    save("data/dqmc_data.jld2",Dict("G"=>G,"τs"=>τs,"β"=>β))
end

A_tmp = zeros(Float64,size(G,2),size(ωs,1))

for k in kpt
    println("DEAC for k ",k)
    out = joinpath(outdir,"$(k).jld2")
    G2 = Matrix{Float64}(G[:,k,:]')
    DEAC_dict = DEAC_Binned(G2,β,τs,ωs,"time_fermionic",100,100,out,"chk.jld2",number_of_generations=20_000,population_size=12,verbose=true)
    A_tmp[k,:] = DEAC_dict["A"][:,end]
end

A_out = zeros(Float64,size(A_tmp,1)+1,size(ωs,1))
A_out[1:16,:] = A_tmp[17:32,:]
A_out[17:33,:] = A_tmp[1:17,:]

out_dict = Dict{String,Any}(
    "A" => A_out,
    "ωs" => ωs
)

save("data/deac.jld2",out_dict)
