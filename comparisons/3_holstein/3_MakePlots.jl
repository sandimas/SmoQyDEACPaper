using CairoMakie
using FileIO
# using ColorShemes

deac = load("data/deac.jld2")
dmrg = load("data/dmrg.jld2")

ks = dmrg["ks"]

f = Figure(size=(2000,800),fontsize=50)
ga = f[1,1] = GridLayout()
maxcol = max(maximum(deac["A"]),maximum(dmrg["A"]) )

ax_dmrg = Axis(ga[1,1],ylabel=L"ω/t",xtickalign=1,xtickcolor=:white,yticksize=20,ytickalign=1,ytickcolor=:white,xticksize=20)
# For some reason the dmrg data given to me was normalized to 0.1
heatmap!(ax_dmrg,ks[1:17],dmrg["ωs"],10.0 .* dmrg["A"][1:17,:],colorrange=(0,maxcol) ,colormap=:terrain)

ax_deac = Axis(ga[1,2],yticklabelsvisible=false,xtickalign=1,xtickcolor=:white,yticksize=20,ytickalign=1,ytickcolor=:white,xticksize=20)#,title="DEAC n=1.00")
heatmap!(ax_deac,ks[17:33],deac["ωs"],deac["A"][1:17,:],colorrange=(0,maxcol) ,colormap=:terrain)
linkyaxes!(ax_dmrg,ax_deac)
colgap!(ga,0)

xtick_pos = [-π,-π/2,0,π/2,π]
xtick_labels = [L"-π/a",L"-π/2a",L"0",L"π/2a",L"π/a"]
ytick_pos = [-4,-2,0,2,4]
ytick_labels = [L"-4",L"-2",L"0",L"2",L"4" ]
ax_deac.xticks = (xtick_pos,xtick_labels)
ax_dmrg.xticks = (xtick_pos,xtick_labels)
ax_dmrg.yticks = (ytick_pos,ytick_labels)
ylims!(ax_deac,(-4.0,4.0) )
hl_pos = 0.41401680839999955 + 1.0
hlines!(ax_deac,[-hl_pos,hl_pos],color=:grey,linewidth=3)
hlines!(ax_dmrg,[-hl_pos,hl_pos],color=:grey,linewidth=3)
text!(ax_dmrg,-π,3.3,text=L"\text{\textbf {(a)} DMRG}",color=:white,fontsize=50)
text!(ax_deac,0.0,3.3,text=L"\text{\textbf{(b)} DEAC}",color=:white,fontsize=50)
Label(ga[1,:,Bottom()], L"k", valign = :bottom,padding=(0.0,0.0,0.0,50.0))
f

save("deac_vs_dmrg.png",f)
