using Interpolations
using JLD2
using DelimitedFiles: readdlm
using StructArrays: StructArray

const beauty_kernel_path = "C:\\Users\\feder\\.julia\\dev\\Fluidum\\src\\kernels\\beauty_reso_dubla_Kj\\"
const light_kernel_path = "C:\\Users\\feder\\.julia\\dev\\Fluidum\\src\\kernels\\FastReso_kernels\\"
struct particle{S,T,U}
    name::S
    mass::T
    degeneracy::U
    thermal_kernel_path::S
    total_kernel_path::S
end

struct kernels{T}
    K1eq::T
    K2eq::T
    K1piphi::T
    K2piphi::T
    K1pieta::T
    K2pieta::T	
    K1bulk::T	
    K2bulk::T
end

struct kernels_diff{T}
    K1eq::T
    K2eq::T
    K1piphi::T
    K2piphi::T
    K1pieta::T
    K2pieta::T	
    K1bulk::T	
    K2bulk::T
    K1diff::T	
    K2diff::T
end

struct kernel_equilibrium{T}
    K1eq::T
    K2eq::T
end

function LightFlavours(;name_file=string("C:\\Users\\feder\\.julia\\dev\\Fluidum\\src\\kernels\\FastReso_kernels/light.txt"),Maxmass=0.939,Minmass=0.0,condition=x->true)
    data     =readdlm(name_file,comment_char='#',comments=true)
    names    =convert.(String,data[:,1])
    mass     =convert.(Float64,data[:,2])
    spin     =convert.(Float64,data[:,4])
    thermal_kernel_path=[] #inizializza un array vuoto 
    [push!(thermal_kernel_path,string(light_kernel_path*names[i]*"_thermal_1565.dat")) for i in eachindex(names)] #tanti path quanti sono i names
    total_kernel_path=[]
    [push!(total_kernel_path,string(light_kernel_path*names[i]*"_total_1565.dat")) for i in eachindex(names)]
    
    return StructArray(particle.(
        names,
        mass,
        2*spin.+1,
        thermal_kernel_path,
        total_kernel_path
       ))

end



function Beauty(;name_file=string("C:\\Users\\feder\\.julia\\dev\\Fluidum\\src\\kernels\\beauty_reso_dubla_Kj\\OpenBeautyParticleList_corrJS.txt"),Maxmass=5.,Minmass=11.,condition=x->true)
    data     =readdlm(name_file,comment_char='#',comments=true)
    names    =convert.(String,data[:,1])
    mass     =convert.(Float64,data[:,2])

    degeneracy=convert.(Float64,data[:,4])

    thermal_kernel_path=[]
    [push!(thermal_kernel_path,string(beauty_kernel_path*names[i]*"_thermal_T0.1565_Kj.out")) for i in eachindex(names)]
    total_kernel_path=[]
    [push!(total_kernel_path,string(beauty_kernel_path*names[i]*"_total_T0.1565_Kj.out")) for i in eachindex(names)]
    
    return StructArray(particle.(
        names,
        mass,
        degeneracy,
        thermal_kernel_path,
        total_kernel_path
       ))
end
""""
returns kernel structure for protons, pions, and kaons. Each object in the structure is an int. func. of (pt,ur)
"""

function save_thermal_kernels()
    
    ext = [] #svuoto l array
    B5279plu,B5324plu,Bs5366zer,Ups9460zer,Ups10000zer = Beauty()

    for i in 1:5 
    part = [B5279plu,B5324plu,Bs5366zer,Ups9460zer,Ups10000zer][i]
    a = readdlm(part.thermal_kernel_path, Float64,skipstart=1)
    #same for all
    pt_list=reshape(a[:,1],(length(unique(a[:,1])),length(unique(a[:,2]))))
    ur_list=reshape(a[:,2],(length(unique(a[:,1])),length(unique(a[:,2]))))
    ur=ur_list[1,:]
    pt=pt_list[:,1]
    itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ]
    push!(ext,extrapolate.(itp, Ref(Flat())) )
    end
    
    B5279plu_thermal=kernels(ext[1]...)
    B5324plu_thermal=kernels(ext[2]...)
    Bs5366zer_thermal=kernels(ext[3]...)
    Ups9460zer_thermal=kernels(ext[4]...)
    Ups10000zer_thermal=kernels(ext[5]...)
    
        #@save "savedfunction.jld" ext
    @save beauty_kernel_path*"B5279plu_thermal.jld" B5279plu_thermal
    @save beauty_kernel_path*"B5324plu_thermal.jld" B5324plu_thermal
    @save beauty_kernel_path*"Bs5366zer_thermal.jld" Bs5366zer_thermal
    @save beauty_kernel_path*"Ups9460zer_thermal.jld" Ups9460zer_thermal
    @save beauty_kernel_path*"Ups10000zer_thermal.jld" Ups10000zer_thermal

    ext = []
    protons,kaons,pions = LightFlavours()
    a = readdlm(protons.thermal_kernel_path, Float64,skipstart=1)
    #same for all
    pt_list=reshape(a[:,1],(length(unique(a[:,1])),length(unique(a[:,2]))))
    ur_list=reshape(a[:,2],(length(unique(a[:,1])),length(unique(a[:,2]))))

    ur=ur_list[1,:]
    pt=pt_list[:,1]

    itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    push!(ext,extrapolate.(itp, Ref(Flat())) )
    a = readdlm(kaons.thermal_kernel_path, Float64,skipstart=1)
    itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    push!(ext,extrapolate.(itp, Ref(Flat())) )
    a = readdlm(pions.thermal_kernel_path, Float64,skipstart=1)
    itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    push!(ext,extrapolate.(itp, Ref(Flat())) )

    protons_thermal=kernels(ext[1]...)
    kaons_thermal=kernels(ext[2]...)
    pions_thermal=kernels(ext[3]...)
    #@save "savedfunction.jld" ext
    @save light_kernel_path*"proton_thermal.jld" protons_thermal
    @save light_kernel_path*"kaon_thermal.jld" kaons_thermal
    @save light_kernel_path*"pion_thermal.jld" pions_thermal
end
#save_thermal_kernels()

function save_total_kernels()
    ext = [] #svuoto l array
    B5279plu,B5324plu,Bs5366zer,Ups9460zer,Ups10000zer = Beauty()

    for i in 1:5 
    part = [B5279plu,B5324plu,Bs5366zer,Ups9460zer,Ups10000zer][i]
    a = readdlm(part.total_kernel_path, Float64,skipstart=1)
    #same for all
    pt_list=reshape(a[:,1],(length(unique(a[:,1])),length(unique(a[:,2]))))
    ur_list=reshape(a[:,2],(length(unique(a[:,1])),length(unique(a[:,2]))))
    ur=ur_list[1,:]
    pt=pt_list[:,1]
    itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ]
    push!(ext,extrapolate.(itp, Ref(Flat())) )
    end
    
    B5279plu_total=kernels(ext[1]...)
    B5324plu_total=kernels(ext[2]...)
    Bs5366zer_total=kernels(ext[3]...)
    Ups9460zer_total=kernels(ext[4]...)
    Ups10000zer_total=kernels(ext[5]...)
    
        #@save "savedfunction.jld" ext
    @save beauty_kernel_path*"B5279plu_total.jld" B5279plu_total
    @save beauty_kernel_path*"B5324plu_total.jld" B5324plu_total
    @save beauty_kernel_path*"Bs5366zer_total.jld" Bs5366zer_total
    @save beauty_kernel_path*"Ups9460zer_total.jld" Ups9460zer_total
    @save beauty_kernel_path*"Ups10000zer_total.jld" Ups10000zer_total
   

    ext = []
    protons,kaons,pions = LightFlavours()
    a = readdlm(protons.total_kernel_path, Float64,skipstart=1)
    #same for all
    pt_list=reshape(a[:,1],(length(unique(a[:,1])),length(unique(a[:,2]))))
    ur_list=reshape(a[:,2],(length(unique(a[:,1])),length(unique(a[:,2]))))

    ur=ur_list[1,:]
    pt=pt_list[:,1]

    itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    push!(ext,extrapolate.(itp, Ref(Flat())) )
    a = readdlm(kaons.total_kernel_path, Float64,skipstart=1)
    itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    push!(ext,extrapolate.(itp, Ref(Flat())) )
    a = readdlm(pions.total_kernel_path, Float64,skipstart=1)
    itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    push!(ext,extrapolate.(itp, Ref(Flat())) )

    protons_total=kernels(ext[1]...)
    kaons_total=kernels(ext[2]...)
    pions_total=kernels(ext[3]...)
    #@save "savedfunction.jld" ext
    @save light_kernel_path*"proton_total.jld" protons_total
    @save light_kernel_path*"kaon_total.jld" kaons_total
    @save light_kernel_path*"pion_total.jld" pions_total
    #return (protons=kernels(ext[1]...),kaons=kernels(ext[2]...),pions=kernels(ext[3]...))
    
    # ext = []
    # D0,Dplus,Dstar0,Dstarplus,Dsplus,Lcplus, Omc, Xic = OpenCharm()
    # a = readdlm(D0.total_kernel_path, Float64,skipstart=1)
    # #same for all
    # pt_list=reshape(a[:,1],(length(unique(a[:,1])),length(unique(a[:,2]))))
    # ur_list=reshape(a[:,2],(length(unique(a[:,1])),length(unique(a[:,2]))))

    # ur=ur_list[1,:]
    # pt=pt_list[:,1]

    # itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    # push!(ext,extrapolate.(itp, Ref(Flat())) )
    # a = readdlm(Dplus.total_kernel_path, Float64,skipstart=1)
    # itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    # push!(ext,extrapolate.(itp, Ref(Flat())) )
    # a = readdlm(Dstar0.total_kernel_path, Float64,skipstart=1)
    # itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    # push!(ext,extrapolate.(itp, Ref(Flat())) )
    # a = readdlm(Dstarplus.total_kernel_path, Float64,skipstart=1)
    # itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    # push!(ext,extrapolate.(itp, Ref(Flat())) )
    # a = readdlm(Dsplus.total_kernel_path, Float64,skipstart=1)
    # itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    # push!(ext,extrapolate.(itp, Ref(Flat())) )
    # a = readdlm(Lcplus.total_kernel_path, Float64,skipstart=1)
    # itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    # push!(ext,extrapolate.(itp, Ref(Flat())) )
    # a = readdlm(Omc.total_kernel_path, Float64,skipstart=1)
    # itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    # push!(ext,extrapolate.(itp, Ref(Flat())) )
    # a = readdlm(Xic.total_kernel_path, Float64,skipstart=1)
    # itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    # push!(ext,extrapolate.(itp, Ref(Flat())) )



    # D0_total=kernels_diff(ext[1]...)
    # Dplus_total=kernels_diff(ext[2]...)
    # Dstar0_total=kernels_diff(ext[3]...)
    # Dstarplus_total=kernels_diff(ext[4]...)
    # Dsplus_total=kernels_diff(ext[5]...)
    # Lcplus_total=kernels_diff(ext[6]...)
    # Omc_total=kernels_diff(ext[7]...)
    # Xic_total=kernels_diff(ext[8]...)
    
    # #@save "savedfunction.jld" ext
    # @save string(@__DIR__,"/Kernels_folder/FastReso_OC_kernels/D0_total.jld") D0_total
    # @save string(@__DIR__,"/Kernels_folder/FastReso_OC_kernels/Dplus_total.jld") Dplus_total
    # @save string(@__DIR__,"/Kernels_folder/FastReso_OC_kernels/Dstar0_total.jld") Dstar0_total
    # @save string(@__DIR__,"/Kernels_folder/FastReso_OC_kernels/Dstarplus_total.jld") Dstarplus_total
    # @save string(@__DIR__,"/Kernels_folder/FastReso_OC_kernels/Dsplus_total.jld") Dsplus_total
    # @save string(@__DIR__,"/Kernels_folder/FastReso_OC_kernels/Lcplus_total.jld") Lcplus_total
    # @save string(@__DIR__,"/Kernels_folder/FastReso_OC_kernels/Omc_total.jld") Omc_total
    # @save string(@__DIR__,"/Kernels_folder/FastReso_OC_kernels/Xic_total.jld") Xic_total


    # ext = []
    # Jpsi = HiddenCharm()[1]
    # a = readdlm(Jpsi.total_kernel_path, Float64,skipstart=1)
    # #same for all
    # pt_list=reshape(a[:,1],(length(unique(a[:,1])),length(unique(a[:,2]))))
    # ur_list=reshape(a[:,2],(length(unique(a[:,1])),length(unique(a[:,2]))))

    # ur=ur_list[1,:]
    # pt=pt_list[:,1]

    # itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    # push!(ext,extrapolate.(itp, Ref(Flat())) )
    # Jpsi_total=kernels_diff(ext[1]...)
    # @save string(@__DIR__,"/Kernels_folder/FastReso_HC_kernels/Jpsi_total.jld") Jpsi_total
end


#uncomment only if need to generate new kernel extrapolations
#save_thermal_kernels()
#save_total_kernels()

