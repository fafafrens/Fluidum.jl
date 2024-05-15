using Interpolations
using JLD2
using DelimitedFiles: readdlm
using StructArrays: StructArray

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

function LightFlavours(;name_file=string(@__DIR__,"/Kernels_only_01445/FastReso_kernels/particles.data"),Maxmass=0.939,Minmass=0.13,condition=x->true)
    data     =readdlm(name_file,comment_char='#',comments=true)
    names    =convert.(String,data[:,1])
    mass     =convert.(Float64,data[:,2])
    spin     =convert.(Float64,data[:,4])
    thermal_kernel_path=[] #inizializza un array vuoto 
    [push!(thermal_kernel_path,string(@__DIR__,"/Kernels_only_01445/FastReso_kernels/"*names[i]*"_thermal_T0.1445_Kj.out")) for i in eachindex(names)] #tanti path quanti sono i names
    total_kernel_path=[]
    [push!(total_kernel_path,string(@__DIR__,"/Kernels_only_01445/FastReso_kernels/"*names[i]*"_total_T0.1445_Kj.out")) for i in eachindex(names)]
    
    return StructArray(particle.(
        names,
        mass,
        2*spin.+1,
        thermal_kernel_path,
        total_kernel_path
       ))

end

function OpenCharm(;name_file=string(@__DIR__,"/Kernels_folder/FastReso_OC_kernels/particles.data"),Maxmass=0.939,Minmass=0.13,condition=x->true)
    data     =readdlm(name_file,comment_char='#',comments=true)
    names    =convert.(String,data[:,1])
    mass     =convert.(Float64,data[:,2])

    degeneracy=convert.(Float64,data[:,4])

    thermal_kernel_path=[]
    [push!(thermal_kernel_path,string(@__DIR__,"/Kernels_folder/FastReso_OC_kernels/"*names[i]*"_thermal_T0.1565_Kj.out")) for i in eachindex(names)]
    total_kernel_path=[]
    [push!(total_kernel_path,string(@__DIR__,"/Kernels_folder/FastReso_OC_kernels/"*names[i]*"_total_T0.1565_Kj.out")) for i in eachindex(names)]
    
    return StructArray(particle.(
        names,
        mass,
        degeneracy,
        thermal_kernel_path,
        total_kernel_path
       ))
end

function HiddenCharm(;name_file=string(@__DIR__,"/Kernels_folder/FastReso_HC_kernels/particles.data"),Maxmass=0.939,Minmass=0.13,condition=x->true)
    data     =readdlm(name_file,comment_char='#',comments=true)
    names    =convert.(String,data[:,1])
    mass     =convert.(Float64,data[:,2])

    degeneracy=convert.(Float64,data[:,4])

    thermal_kernel_path=[]
    [push!(thermal_kernel_path,string(@__DIR__,"/Kernels_folder/FastReso_HC_kernels/jp3096zer_thermal_T0.1565_Kj.out")) for i in eachindex(names)]
    total_kernel_path=[]
    [push!(total_kernel_path,string(@__DIR__,"/Kernels_folder/FastReso_HC_kernels/jp3096zer_total_T0.1565_Kj.out")) for i in eachindex(names)]
    
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
    ext = []
    protons,kaons,pions = LightFlavours()
    a = readdlm(protons.thermal_kernel_path, Float64,skipstart=1)
    #same for all
    pt_list=reshape(a[:,1],(length(unique(a[:,1])),length(unique(a[:,2])))) #riformatto in modo che pt sia una matrice con colonne tutte uguali e righe con pt 
    ur_list=reshape(a[:,2],(length(unique(a[:,1])),length(unique(a[:,2])))) #riformatto in modo che ur sia una matrice con righe tutte uguali e righe con ur 

    ur=ur_list[1,:]
    pt=pt_list[:,1]

    itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] #interpolazione di tutte le colonne da 3 in poi, tramite iterpolazione lineare
    push!(ext,extrapolate.(itp, Ref(Flat())) ) #primo elemento pushato in ext: sara ext[1]. Estrapolazione con una valore costante dappertutto (Ref(Flat)) Per fare l estrapolazione ci serve una funzione quindi prima ho dovuto interpolare
    a = readdlm(kaons.thermal_kernel_path, Float64,skipstart=1)
    itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    push!(ext,extrapolate.(itp, Ref(Flat())) ) #secondo elemento pushato in ext: sara ext[1]
    a = readdlm(pions.thermal_kernel_path, Float64,skipstart=1)
    itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    push!(ext,extrapolate.(itp, Ref(Flat())) )

    protons_thermal=kernels(ext[1]...) #... indica che sto eseguendo un "unpaking": sto prendendo tutti gli elementi di ext[1], e li sto usando per riempire kernels 
    print("the type: ", typeof(ext[1]))
    #print("k1 eq: ",K1eq.kernels)
    #print("k1 eq: ", K1eq(pt,ur))
    kaons_thermal=kernels(ext[2]...)
    pions_thermal=kernels(ext[3]...)
    #@save "savedfunction.jld" ext
    @save string(@__DIR__,"/Kernels_only_01445/FastReso_kernels/proton_thermal.jld") protons_thermal
    @save string(@__DIR__,"/Kernels_only_01445/FastReso_kernels/kaon_thermal.jld") kaons_thermal
    @save string(@__DIR__,"/Kernels_only_01445/FastReso_kernels/pion_thermal.jld") pions_thermal
    #return (protons=kernels(ext[1]...),kaons=kernels(ext[2]...),pions=kernels(ext[3]...))
    
    
    # ext = [] #svuoto l array
    # D0,Dplus,Dstar0,Dstarplus,Dsplus,Lcplus,Omc,Xic = OpenCharm()
    # a = readdlm(D0.thermal_kernel_path, Float64,skipstart=1)
    # #same for all
    # pt_list=reshape(a[:,1],(length(unique(a[:,1])),length(unique(a[:,2]))))
    # ur_list=reshape(a[:,2],(length(unique(a[:,1])),length(unique(a[:,2]))))

    # ur=ur_list[1,:]
    # pt=pt_list[:,1]

    # itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    # push!(ext,extrapolate.(itp, Ref(Flat())) )
    # a = readdlm(Dplus.thermal_kernel_path, Float64,skipstart=1)
    # itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    # push!(ext,extrapolate.(itp, Ref(Flat())) )
    # a = readdlm(Dstar0.thermal_kernel_path, Float64,skipstart=1)
    # itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    # push!(ext,extrapolate.(itp, Ref(Flat())) )
    # a = readdlm(Dstarplus.thermal_kernel_path, Float64,skipstart=1)
    # itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    # push!(ext,extrapolate.(itp, Ref(Flat())) )
    # a = readdlm(Dsplus.thermal_kernel_path, Float64,skipstart=1)
    # itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    # push!(ext,extrapolate.(itp, Ref(Flat())) )
    # a = readdlm(Lcplus.thermal_kernel_path, Float64,skipstart=1)
    # itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    # push!(ext,extrapolate.(itp, Ref(Flat())) )
    # a = readdlm(Omc.thermal_kernel_path, Float64,skipstart=1)
    # itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    # push!(ext,extrapolate.(itp, Ref(Flat())) )
    # a = readdlm(Xic.thermal_kernel_path, Float64,skipstart=1)
    # itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    # push!(ext,extrapolate.(itp, Ref(Flat())) )


    # D0_thermal=kernels_diff(ext[1]...)
    # Dplus_thermal=kernels_diff(ext[2]...)
    # Dstar0_thermal=kernels_diff(ext[3]...)
    # Dstarplus_thermal=kernels_diff(ext[4]...)
    # Dsplus_thermal=kernels_diff(ext[5]...)
    # Lcplus_thermal=kernels_diff(ext[6]...)
    # Omc_thermal=kernels_diff(ext[7]...)
    # Xic_thermal=kernels_diff(ext[8]...)
    
    # #@save "savedfunction.jld" ext
    # @save string(@__DIR__,"/Kernels_folder/FastReso_OC_kernels/D0_thermal.jld") D0_thermal
    # @save string(@__DIR__,"/Kernels_folder/FastReso_OC_kernels/Dplus_thermal.jld") Dplus_thermal
    # @save string(@__DIR__,"/Kernels_folder/FastReso_OC_kernels/Dstar0_thermal.jld") Dstar0_thermal
    # @save string(@__DIR__,"/Kernels_folder/FastReso_OC_kernels/Dstarplus_thermal.jld") Dstarplus_thermal
    # @save string(@__DIR__,"/Kernels_folder/FastReso_OC_kernels/Dsplus_thermal.jld") Dsplus_thermal
    # @save string(@__DIR__,"/Kernels_folder/FastReso_OC_kernels/Lcplus_thermal.jld") Lcplus_thermal
    # @save string(@__DIR__,"/Kernels_folder/FastReso_OC_kernels/Omc_thermal.jld") Omc_thermal
    # @save string(@__DIR__,"/Kernels_folder/FastReso_OC_kernels/Xic_thermal.jld") Xic_thermal
    

    # ext = []
    # Jpsi = HiddenCharm()[1]
    # a = readdlm(Jpsi.thermal_kernel_path, Float64,skipstart=1)
    # #same for all
    # pt_list=reshape(a[:,1],(length(unique(a[:,1])),length(unique(a[:,2]))))
    # ur_list=reshape(a[:,2],(length(unique(a[:,1])),length(unique(a[:,2]))))

    # ur=ur_list[1,:]
    # pt=pt_list[:,1]

    # itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    # itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ] 
    # push!(ext,extrapolate.(itp, Ref(Flat())) )
    
    # Jpsi_thermal=kernels_diff(ext[1]...)
    # @save string(@__DIR__,"/Kernels_folder/FastReso_HC_kernels/Jpsi_thermal.jld") Jpsi_thermal
end
#save_thermal_kernels()

function save_total_kernels()
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
    @save string(@__DIR__,"/Kernels_only_01445/FastReso_kernels/proton_total.jld") protons_total
    @save string(@__DIR__,"/Kernels_only_01445/FastReso_kernels/kaon_total.jld") kaons_total
    @save string(@__DIR__,"/Kernels_only_01445/FastReso_kernels/pion_total.jld") pions_total
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

