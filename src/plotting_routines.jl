#=function plot_params(;gui = false,pHeight=2.3)
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
=#
using DelimitedFiles

# Function to process percentage strings
function process_value(value)
    if occursin("%", value)
        # Remove the '%' and divide the number by 100
        return parse(Float64, replace(value, "%"=>"")) / 100
    else
        # Otherwise, return the value as a string or number
        try
            return parse(Float64, value)  # Try to convert to a number if possible
        catch
            return value  # If not, return as a string
        end
    end
end

# Function to read CSV, filter comments, and process percentages
function read_csv_with_header(file_path)
    open(file_path, "r") do file
        lines = filter(line -> !startswith(line, "#"), readlines(file))
        
        # Extract the header (assuming the first non-commented line is the header)
        header = split(lines[1], ",")  # Split the header by commas
        data_lines = lines[2:end]  # Remaining lines are data
        
        # Process each data line
        processed_data_lines = []
        for line in data_lines
            columns = split(line, ",")
            processed_columns = [process_value(col) for col in columns]  # Process each value
            push!(processed_data_lines, join(processed_columns, ","))
        end
        
        # Write the processed data (without # lines) to a temporary file
        open("filtered_file.csv", "w") do out_file
            write(out_file, join(processed_data_lines, "\n"))
        end

        # Read the filtered and processed data
        data = readdlm("filtered_file.csv", ',')
        
        return header, data
    end
end

function plot_spectra(obj::Observables{S,T,U,V,M,K,A,B,C,D};path="./plots/",thermal=true,total=false,save=false,norm_spectra=1.) where {S,T,U,V,M,K,A,B,C,D}
    if isdir(path)
    else mkdir(path)
    end
    #fig1, ax1 = subplots()
    name = obj.particle.name
    title = L"\mathrm{non-prompt}\, \mathrm{D_0}"
    if thermal == true
        plot(obj.pt_bins,2π*norm_spectra *obj.pt_bins .*obj.spectra_th,c="#FFC20A",label=L"\mathrm{w/o\, res.\, dec.}",yaxis=:log)
    end
    if total == true
        plot(obj.pt_bins,2π*norm_spectra *obj.pt_bins .*obj.spectra_tot,c="#FFC20A",label=L"\mathrm{w/\, res.\, dec.}",yaxis=:log,xlims=(0,10),legend_title=title, frame=true, grid=false, size=(250*1.618,250),legend_font_pointsize=8,legend_title_font_pointsize=8)
    end
    if name == "D0"
        title = L"\mathrm{non-prompt}\, \mathrm{D_0}"
        header, data = read_csv_with_header("C:\\Users\\feder\\Downloads\\HEPData-ins2025044-v1-Table_1a.csv")
    elseif name == "D_s+"
        title = L"\mathrm{non-prompt}\, \mathrm{D_s^+}"
        header, data = read_csv_with_header("C:\\Users\\feder\\Downloads\\HEPData-ins2071181-v1-Table_0.csv")
    elseif name == "Jpsi"
        title = L"\mathrm{non-prompt}\, \mathrm{J/\psi}"
        header, data = read_csv_with_header("C:\\Users\\feder\\Downloads\\HEPData-ins2692201-v1-Table_5.csv")
    end
    pt_bin_central = data[:,1]
    pt_bin_low = data[:,2]
    pt_bin_high = data[:,3]
    dndpt = data[:,4]
    err_plu = zeros(length(dndpt))
    err_min = zeros(length(dndpt))
    if length(dndpt)>8
        [data[:,i]*=dndpt[i] for i in 9:lastindex(data[1,:])]
    end
    [err_plu .+=data[:,i].^2 for i in range(5,lastindex(data[1,:]),step=2)]
    [err_min .+=data[:,i].^2 for i in range(6,lastindex(data[1,:]),step=2)]
    err_plu .= sqrt.(err_plu)
    err_min .= sqrt.(err_min)
    #err_plu = sqrt.(data[:,5].^2 .+data[:,7].^2)
    #err_min =sqrt.(data[:,6].^2 .+data[:,8].^2)
    yerr = [(err_min[i],err_plu[i]) for i in eachindex(err_min)]
    xerr = [(pt_bin_high[i]-pt_bin_low[i])*0.5 for i in eachindex(pt_bin_low)]
    plot!(pt_bin_central,dndpt,yerr=yerr,xerr=xerr,label=L"\mathrm{ALICE}",c="black")
    xlabel!(L"\mathrm{p_T\,[GeV]}")
    ylabel!(L"\mathrm{dN/dp_Tdy\,[GeV^{-1}]}")

    #yticks!([10e-6,10e-5,10e-4,10e-3,10e-2,10e-1,1])
    #fig1.set_size_inches(pWidth, pHeight)
    if save==true
        filename = Fluidum.get_filename(obj,path=path,toplot=true)
        #fig1.savefig(filename,bbox_inches="tight")
        savefig(filename)
    end
end

function plot_spectra(obj::Observables{S,T,U,V,M,K,A,B,C,D},obj2::Observables{S,T,U,V,M,K,A,B,C,D};path="./plots/",thermal=true,total=false,save=false,norm_spectra1=1.,norm_spectra2=1.) where {S,T,U,V,M,K,A,B,C,D}
    if isdir(path)
    else mkdir(path)
    end
    name = obj.particle.name
    title = L"\mathrm{"*name*L"}"
        
    #fig1, ax1 = subplots()
    if thermal == true
        mid = 2π*obj.pt_bins .*(norm_spectra1 *obj.spectra_th .+ norm_spectra2 *obj2.spectra_th) ./ 2   #the midpoints (usually representing mean values)
        w = 2π*norm_spectra1 *obj.pt_bins .*(norm_spectra2 *obj2.spectra_th .- norm_spectra1 *obj.spectra_th) ./ 2  
        plot(obj.pt_bins, mid, ribbon = w,c="#FFC20A",label=L"\mathrm{w/o\, res.\, dec.}",yaxis=:log,xlims=(0,10),legend_title=title, framestyle = :box, grid=false, size=(250*1.618,250),legend_font_pointsize=8,legend_title_font_pointsize=8,legend=:bottomright) 
    
        #plot(obj.pt_bins,2π*norm_spectra1 *obj.pt_bins .*obj.spectra_th,fillrange=2π*norm_spectra2 *obj.pt_bins .*obj2.spectra_th,fillalpha = 0.5,c="#FFC20A",label=L"\mathrm{w/o\, res.\, dec.}",yaxis=:log)
    end
    if total == true
        mid = 2π*obj.pt_bins .*(norm_spectra1 *obj.spectra_tot .+ norm_spectra2 *obj2.spectra_tot) ./ 2   #the midpoints (usually representing mean values)
        w = 2π*norm_spectra1 *obj.pt_bins .*(norm_spectra2 *obj2.spectra_tot .- norm_spectra1 *obj.spectra_tot) ./ 2  
        plot!(obj.pt_bins, mid, ribbon = w,c="#0C7BDC",label=L"\mathrm{w/\, res.\, dec.}",yaxis=:log,xlims=(0,10),legend_title=title, framestyle = :box, grid=false, size=(250*1.618,250),legend_font_pointsize=8,legend_title_font_pointsize=8,legend=:bottomright) 
    
       # plot(obj.pt_bins,2π*norm_spectra1 *obj.pt_bins .*obj.spectra_tot,fillrange=2π*norm_spectra2 *obj.pt_bins .*obj2.spectra_tot,fillalpha = 0.5,c="#0C7BDC",label=L"\mathrm{w/\, res.\, dec.}",yaxis=:log,xlims=(0,10),legend_title=title, framestyle = :box, grid=false, size=(250*1.618,250),legend_font_pointsize=8,legend_title_font_pointsize=8,legend=:bottomright)
    end
    name = obj.particle.name
    exp = false
    if name == "D0"
        exp = true
        title = L"\mathrm{non-prompt}\, \mathrm{D_0}"
        header, data = read_csv_with_header("C:\\Users\\feder\\Downloads\\HEPData-ins2025044-v1-Table_1a.csv")
    elseif name == "D_s+"
        exp = true
        title = L"\mathrm{non-prompt}\, \mathrm{D_s^+}"
        header, data = read_csv_with_header("C:\\Users\\feder\\Downloads\\HEPData-ins2071181-v1-Table_0.csv")
    elseif name == "Jpsi"
        exp = true
        title = L"\mathrm{non-prompt}\, \mathrm{J/\psi}"
        header, data = read_csv_with_header("C:\\Users\\feder\\Downloads\\HEPData-ins2692201-v1-Table_5.csv")
    end

    if exp == true
    pt_bin_central = data[:,1]
    pt_bin_low = data[:,2]
    pt_bin_high = data[:,3]
    dndpt = data[:,4]
    err_plu = zeros(length(dndpt))
    err_min = zeros(length(dndpt))
    if length(dndpt)>8
        [data[:,i]*=dndpt[i] for i in 9:lastindex(data[1,:])]
    end
    #[err_plu .+=data[:,i].^2 for i in range(5,lastindex(data[1,:]),step=2)]
    #[err_min .+=data[:,i].^2 for i in range(6,lastindex(data[1,:]),step=2)]
    #err_plu .= sqrt.(err_plu)
    #err_min .= sqrt.(err_min)
    err_plu = sqrt.(data[:,5].^2 .+data[:,7].^2)
    err_min =sqrt.(data[:,6].^2 .+data[:,8].^2)
    yerr = [(err_min[i],err_plu[i]) for i in eachindex(err_min)]
    xerr = [(pt_bin_high[i]-pt_bin_low[i])*0.5 for i in eachindex(pt_bin_low)]
    plot!(pt_bin_central,dndpt,yerr=yerr,xerr=xerr,label=L"\mathrm{ALICE}",c="black",xlims=(0,10),legend_title=title, grid=false, size=(250*1.618,250),legend_font_pointsize=8,legend_title_font_pointsize=8)
    end
    xlabel!(L"\mathrm{p_T\,[GeV]}")
    ylabel!(L"\mathrm{dN/dp_Tdy\,[GeV^{-1}]}")

    #yticks!([10e-6,10e-5,10e-4,10e-3,10e-2,10e-1,1])
    #fig1.set_size_inches(pWidth, pHeight)
    if save==true
        filename = Fluidum.get_filename(obj,path=path,toplot=true)
        #fig1.savefig(filename,bbox_inches="tight")
        savefig(filename)
    end
end

function plot_int_yields(obs1::Observables{S,T,U,V,M,K,A,B,C,D},obs2::Observables{S,T,U,V,M,K,A,B,C,D})  where {S,T,U,V,M,K,A,B,C,D}
    x = obs1.particle.mass
    y1 = obs1.yield_th
    y2 = obs2.yield_th
    mid = (y1 .+ y2)/2
    w = (y2 .- y1)/2
    
    plot(x,mid, ribbon = w,c="#FFC20A",label=L"\mathrm{w/o\, res.\, dec.}")

    y1 = obs1.yield_tot
    y2 = obs2.yield_tot
    mid = (y1 .+ y2)/2
    w = (y2 .- y1)/2
    plot!(x,mid, ribbon = w,c="#0C7BDC",label=L"\mathrm{w/\, res.\, dec.}")
end

function plot_int_yields!(obs1::Observables{S,T,U,V,M,K,A,B,C,D},obs2::Observables{S,T,U,V,M,K,A,B,C,D})  where {S,T,U,V,M,K,A,B,C,D}
    x = obs1.particle.mass
    y1 = obs1.yield_th
    y2 = obs2.yield_th
    mid = (y1 .+ y2)/2
    w = (y2 .- y1)/2
    
    plot!(x,mid, ribbon = w,c="#FFC20A",label=L"\mathrm{w/o\, res.\, dec.}")

    y1 = obs1.yield_tot
    y2 = obs2.yield_tot
    mid = (y1 .+ y2)/2
    w = (y2 .- y1)/2
    plot!(x,mid, ribbon = w,c="#0C7BDC",label=L"\mathrm{w/\, res.\, dec.}")
end

function plot_int_yields(obs1,obs2;save=true)
    x = zeros(length(obs1))
    y1 = zeros(length(obs1))
    y2 = zeros(length(obs1))
    
    for i in eachindex(obs1)
    x[i] = obs1[i].particle.mass
    y1[i] = obs1[i].yield_th
    y2[i] = obs2[i].yield_th
    end
    mid = (y1 .+ y2)/2
    w = (y2 .- y1)/2
  
    plot(x,mid, yerror = w,c="#FFC20A",label=L"\mathrm{w/o\, res.\, dec.}",seriestype=:scatter)

    for i in eachindex(obs1)
        x[i] = obs1[i].particle.mass
        y1[i] = obs1[i].yield_tot
        y2[i] = obs2[i].yield_tot
        end
    mid = (y1 .+ y2)/2
    w = (y2 .- y1)/2
    plot!(x,mid, yerror = w,c="#0C7BDC",label=L"\mathrm{w/\, res.\, dec.}",seriestype=:scatter, grid=false,  framestyle = :box, size=(300,250))

    xlabel!(L"\mathrm{M\,[GeV]}")
    ylabel!(L"\mathrm{dN/dy}")

    if save==true
        savefig("yields.pdf")
    end
end

"""
    plot_field(fields,express::Symbol; tspan=(1,10),Δt=2,gridpoints=500,rmax=30)

TBW
"""
function plot_field(fields,express::Symbol; tspan=(1,10),Δt=2,gridpoints=500,rmax=30,save=false)
    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    oned_visc_hydro = Fluidum.HQ_viscous_1d()
    disc_fields = DiscreteFileds(oned_visc_hydro,disc,Float64) 
    index=get_index(express,disc_fields.fields)
    plot_field(fields,index; tspan=tspan,Δt=Δt,gridpoints=gridpoints,rmax=rmax)
    if save==true
        savefig("./plots/$express._1.pdf")
    end
end

function plot_field(fields,index::Int64; tspan=(1,10),Δt=2,gridpoints=500,rmax=30)
    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    x = [disc.grid[i][1] for i in 2:lastindex(disc.grid)-1]
    t = tspan[1]
    plot(x,fields(tspan[1])[index,2:lastindex(disc.grid)-1], xlabel="r [fm]",ylabel="field", label="t = $t [fm]")
    [plot!(x,fields(tspan[1]+i*Δt)[index,2:lastindex(disc.grid)-1], label="t = $(tspan[1]+i*Δt) [fm]") for i in 1:floor((tspan[2]-tspan[1])/Δt)]
end

function plot_field(fields; tspan=(1,10),Δt=2,gridpoints=500,rmax=30,save=false)
    [plot_field(fields,index; tspan=tspan,Δt=Δt,gridpoints=gridpoints,rmax=rmax,save=save) for index in 1:1:length(fields)]
end

function plot_density(result; tspan=(1,10),Δt=2,gridpoints=500,rmax=30,eos=HadronResonaceGas(), save = false)
    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    x = [disc.grid[i][1] for i in 2:lastindex(disc.grid)-1]
    t = tspan[1]
    n(t) = getindex.(getproperty.(thermodynamic.(result(t)[1,2:lastindex(disc.grid)-1],result(t)[6,2:lastindex(disc.grid)-1],Ref(eos)),:pressure),1)
    dndm(t) = getindex.(getproperty.(thermodynamic.(result(t)[1,2:lastindex(disc.grid)-1],result(t)[6,2:lastindex(disc.grid)-1],Ref(eos)),:pressure_derivative),2)
    dndt(t) = getindex.(getproperty.(thermodynamic.(result(t)[1,2:lastindex(disc.grid)-1],result(t)[6,2:lastindex(disc.grid)-1],Ref(eos)),:pressure_derivative),1)

    plot(x,n(tspan[1]), xlabel="r [fm]",ylabel="n", label="t = $t [fm]")
    [plot!(x,n(tspan[1]+i*Δt), label="t = $(tspan[1]+i*Δt) [fm]") for i in 1:floor((tspan[2]-tspan[1])/Δt)]
    if save==true
        savefig("./plots/density.pdf")
    end
end

