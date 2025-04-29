function plot_params(;gui = false,pHeight=2.3)
    close("all")
    pygui(gui)
    rc("text", usetex ="true")
    rc("font", family = "serif")
    rc("font", size = 12)
    rc("xtick", labelsize=12)
    rc("ytick", labelsize=12)
    rc("lines", markersize = 1)
    rc("legend", fontsize = 10)
    rc("axes", titlesize = 12)
    rc("axes", titlepad = 10)
    
    #x = [discretization.grid[i][1] for i in eachindex(x)]
    #rc("keymap.quit")
    #fig, ax = subplots()#figsize=(10,15));
    #fig.set_size_inches(pWidth, pHeight)
    #ax.tick_params(axis="both", which="major")
    return nothing
end 
pHeight = 2.3
pWidth = 1.3333333 * pHeight
#pHeight = 1.3333333 *5
#pWidth = 5
export pHeight,pWidth    

function plot_spectra(obj::Observables{S,T,U,V,M,K,A,B,C,D};path="./plots/",thermal=true,total=false,save=false) where {S,T,U,V,M,K,A,B,C,D}
    if isdir(path)
    else mkdir(path)
    end
    fig1, ax1 = subplots()
    if thermal == true
        ax1.plot(obj.pt_bins,2π *obj.pt_bins .*obj.spectra_th,"#FFC20A",label=L"\mathrm{w/o\, res.\, dec.}")
    end
    if total == true
        ax1.plot(obj.pt_bins,2π *obj.pt_bins .*obj.spectra_tot,"#FFC20A",label=L"\mathrm{w/\, res.\, dec.}")
    end
    
    fig1.set_size_inches(pWidth, pHeight)
    if save==true
        filename = Fluidum.get_filename(obj,path=path,toplot=true)
        fig1.savefig(filename,bbox_inches="tight")
    end
end