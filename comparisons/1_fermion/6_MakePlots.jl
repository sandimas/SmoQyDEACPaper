using CairoMakie
using FileIO
using Makie
using Printf
using LaTeXStrings

folders = ["01_Gauss_1.0","02_Gauss_0.1","03_Gauss_0.01",
            "04_Lorentz_1.0","05_Double_Gauss_1","06_Double_Gauss_2",
            "07_Gauss_Singular","08_Log_Singular","09_Chaos","10_Chaos_Singular"]


function plots(folder)
    name = folder
    
    f= Figure(resolution=(1000,750),fontsize=30)
    g = f[1,1] = GridLayout()
    
    ax1 = Axis(g[1,1],xtickalign=1,xticklabelsvisible=false,ytickalign=1)
    ax2 = Axis(g[1,2],xtickalign=1,xticklabelsvisible=false,yticklabelsvisible=false,ytickalign=1)

    ax3 = Axis(g[2,1],xtickalign=1,ytickalign=1)
    ax4 = Axis(g[2,2],xtickalign=1,ytickalign=1,yticklabelsvisible=false)

    colgap!(g,0)
    rowgap!(g,0)
    linkxaxes!(ax1,ax3)
    linkxaxes!(ax2,ax4)
    linkyaxes!(ax1,ax2)
    linkyaxes!(ax3,ax4)
    try
        dict = load(joinpath(folder,"true.jld2"))
        lines!(ax1,dict["ωs"],dict["A"],color=:black,label="True",linewidth=5,grid=false)
        lines!(ax2,dict["ωs"],dict["A"],color=:black,linewidth=5)
        lines!(ax3,dict["ωs"],dict["A"],color=:black,linewidth=5)
        lines!(ax4,dict["ωs"],dict["A"],color=:black,linewidth=5)
    catch
    end
    try
        dict = load(joinpath(folder,"deac.jld2"))
        band!(ax1,dict["ωs"],dict["A"] .- dict["σ"],dict["A"] .+ dict["σ"],color=(:blue,0.2))
        lines!(ax1,dict["ωs"],dict["A"],color=:blue,label="DEAC",linewidth=3,grid=false)
        dict = load(joinpath(folder,"deac_1.jld2"))
        lines!(ax1,dict["ωs"],dict["A"],color=:orange,label="DEAC",linewidth=3,grid=false)
        
    catch
    end
    try
        dict = load(joinpath(folder,"mem.jld2"))
        lines!(ax3,dict["ωs"],dict["A"] ,color=:red,label="MaxEnt",linewidth=3)
    catch
    end

    try
        dict = load(joinpath(folder,"mem_cov.jld2"))
        lines!(ax4,dict["ωs"],dict["A"] ,color=:sienna,label="MaxEnt Cov",linewidth=3)
    catch
    end

    try
        dict = load(joinpath(folder,"sac.jld2"))
        lines!(ax2,dict["ωs"],dict["A"] ,color=:green,label="SAC",linewidth=3)
    catch
    end
    xs = (0.,3.5)
    xlims!(ax1,xs)
    xlims!(ax2,xs)
    xlims!(ax3,xs)
    xlims!(ax4,xs)

    

    text!(ax1,-0.25,0.75,text="DEAC",fontsize=30)
    text!(ax2,-0.25,0.75,text="SAC",fontsize=30)
    text!(ax3,-0.25,0.75,text="MEM",fontsize=30)
    text!(ax4,-0.25,0.75,text="MEM Cov",fontsize=30)
    Label(g[2,:,Bottom()], "ω/t", valign = :bottom,padding=(0.0,0.0,0.0,30.0))


    Label(g[:,0], "s(ω)",rotation=π/2)
    save(name*".png",f)


end

println("Plots for paper have had their axes shifted to be easier to see than this function will give you.")

for folder in folders
    plot(folder)
end