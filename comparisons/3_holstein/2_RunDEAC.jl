using FileIO
using SmoQyDEAC
include("../9_scripts/SQMCloader.jl")

directory = "new_1D_w1.00_l0.50_n1.00_L32_b40.00-1/"
if Threads.nthreads() == 1
    println("This should be run \$ julia --threads=auto 2_RunDEAC.jl")
end

kpt = 1:32
beta = 40.0



ωmin = -10.0
ωmax = 10.0
nω = 200
ωs = collect(LinRange(ωmin,ωmax,nω))

outdir = "data/DEAC/"
try
    mkdir(outdir)
catch
end
dict,real,imag,_,_,β = load_from_SmoQyDQMC(simulationfolder=directory,
                                         correlation="greens",
                                         space="momentum",
                                         type="time_displaced",bin=true)

G = real[1,1,:,:,1,1,:,1]
τs=collect(LinRange(0.0,β,size(G,1)))

A_tmp = zeros(Float64,size(G,2),size(ωs,1))

for k in kpt
    println("DEAC for k ",k)
    out = joinpath(outdir,string(kpt)*".jld2")
    G2 = G[:,k,:]'
    DEAC_dict = DEAC_Binned(G2,β,τs,ωs,"time_fermionic",100,10,out,directory,number_of_generations=1000000,population_size=12,verbose=true)
    A_tmp[k,:] = DEAC_dict["A"]
end

A_out = zeros(Float64,size(A_tmp,1)+1,size(ωs,1))
A_out[1:16,:] = A_tmp[17:32,:]
A_out[17:33,:] = A_tmp[1:17,:]

out_dict = Dict{String,Any}(
    "A" => A_out,
    "ωs" => ωs
)

save("data/deac.jld2",out_dict)
