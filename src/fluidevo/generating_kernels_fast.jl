using Interpolations
using JLD2
using DelimitedFiles: readdlm
using StructArrays: StructArray

fluidum_dir="C:/Users/feder/.julia/dev/Fluidum/"

non_prompt_kernels = fluidum_dir*"src\\kernels\\non_prompt\\"

struct particle_attribute{S,R,U,V}
    name::S
    mass::R
    degeneracy::U
    charge::U
    thermal_kernel_path::V
    total_kernel_path::V
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

kernels(file) = kernels(file...)

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

kernels_diff(file) = kernels_diff(file...)
struct kernel_equilibrium{T}
    K1eq::T
    K2eq::T
end

function ParticleData(;name_file=non_prompt_kernels*"particles_pythia8_selection.txt",output_path=non_prompt_kernels,kernel_type="non_prompt",format = "jld2",Maxmass=4.,Minmass=0.0,condition=x->true,loading=true)
    data     =readdlm(name_file,comment_char='#',comments=true)
    names    =convert.(String,data[:,1])
    mass     =convert.(Float64,data[:,2])
    deg      =convert.(Float64,data[:,4])
    spin     =convert.(Float64,data[:,5])
    charge   =convert.(Float64,data[:,12]) +convert.(Float64,data[:,13])
    
    thermal_kernel_path=[] #inizializza un array vuoto 
    [push!(thermal_kernel_path,
    if loading 
        output_path*string(names[i])*"_thermal_T0.1565_Kj.out"
    else 
        load_object(output_path*string(names[i])*"_thermal_T0.1565."*format)
    end
    ) for i in eachindex(names)] #tanti path quanti sono i names
    
    total_kernel_path=[]
    [push!(total_kernel_path,
    if loading
        output_path*string(names[i])*"_"*kernel_type*"_T0.1565_Kj.out"
    else
        load_object(output_path*string(names[i])*"_"*kernel_type*"_T0.1565."*format)
    end
    ) for i in eachindex(names)]
    
    fulllist = StructArray(particle_attribute.(
        names,
        mass,
        deg,
        charge,
        thermal_kernel_path,
        total_kernel_path
       ))
    return filter(x->(x.mass<Maxmass&&x.mass>Minmass&&condition(x)),fulllist)

end

"""
returns kernel structure for protons, pions, and kaons. Each object in the structure is an int. func. of (pt,ur)
"""
function save_kernels(;name_file=non_prompt_kernels*"particles_pythia8_selection.txt",output_path=non_prompt_kernels,kernel_type="nonprompt",format="jld2",Maxmass=4.,Minmass=0.0,condition=x->true)
    print("saving kernels...")
    selection=ParticleData(name_file=name_file,output_path=output_path,kernel_type=kernel_type,Maxmass=Maxmass,Minmass=Minmass,format=format,condition=condition,loading=true)
    print(eachindex(selection))
        for i in eachindex(selection)

          
        a = readdlm(selection[i].thermal_kernel_path, Float64,skipstart=1)
        #same for all
        pt_list=reshape(a[:,1],(length(unique(a[:,1])),length(unique(a[:,2]))))
        ur_list=reshape(a[:,2],(length(unique(a[:,1])),length(unique(a[:,2]))))
        ur=ur_list[1,:]
        pt=pt_list[:,1]
        itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ]
        ext=extrapolate.(itp, Ref(Flat())) 
        save_object(output_path*string(selection[i].name)*"_thermal_T0.1565."*format,ext)
print("saved "*output_path*string(selection[i].name)*"_thermal_T0.1565."*format)
        a = readdlm(selection[i].total_kernel_path, Float64,skipstart=1)
        #same for all
        pt_list=reshape(a[:,1],(length(unique(a[:,1])),length(unique(a[:,2]))))
        ur_list=reshape(a[:,2],(length(unique(a[:,1])),length(unique(a[:,2]))))
        ur=ur_list[1,:]
        pt=pt_list[:,1]
        itp=[interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear()) ) for n in 3:size(a)[2] ]
        ext=extrapolate.(itp, Ref(Flat())) 
        save_object(output_path*string(selection[i].name)*"_"*kernel_type*"_T0.1565."*format,ext)
    end
end
##uncomment only if need to generate new kernel extrapolations
#save_kernels(name_file=fluidum_dir*"src/kernels/beauty_reso_dubla_Kj/OpenBeautyParticleList_corrJS.txt",output_path=fluidum_dir*"src/kernels/beauty_reso_dubla_Kj/",format="jld",kernel_type="total",Maxmass=10.1)
#save_kernels(name_file=fluidum_dir*"src/kernels/prompt_beauty/particles_pythia8_selection.txt",output_path=fluidum_dir*"src/kernels/prompt_beauty/",format="jld2",kernel_type="total",Maxmass=10.1)
#save_kernels(name_file=fluidum_dir*"src/kernels/charm_reso_dubla_Kj/particles.data",output_path=fluidum_dir*"src/kernels/charm_reso_dubla_Kj/",format="jld2",kernel_type="total",Maxmass=10.1)

    
function dictionary(;name_file=fluidum_dir*"src/kernels/non_prompt/particles_pythia8_selection.txt",output_path=non_prompt_kernels,format="jld2",kernel_type="nonprompt",Maxmass=10.1,Minmass=0.0,condition=x->true)
    selection=ParticleData(name_file=name_file,kernel_type=kernel_type,Maxmass=Maxmass,Minmass=Minmass,condition=condition,output_path=output_path,format=format, loading = false)
    
    Key_Tuple=[Meta.parse(replace(selection[i].name, "+" => "plu","(" => "_",")" => "_")) for i in eachindex(selection)]
    return (; zip(Key_Tuple, selection)...) 
    
end 

dic_NP = dictionary();
export dic_NP

dic_P = dictionary(name_file=fluidum_dir*"src/kernels/beauty_reso_dubla_Kj/OpenBeautyParticleList_corrJS.txt",output_path=fluidum_dir*"src/kernels/beauty_reso_dubla_Kj/",format="jld",kernel_type="total")
export dic_P

dic_P_pythia = dictionary(name_file=fluidum_dir*"src/kernels/prompt_beauty/particles_pythia8_selection.txt",output_path=fluidum_dir*"src/kernels/prompt_beauty/",format="jld2",kernel_type="total",Maxmass=10.1)
keys(dic_P_pythia)

dic_charm = dictionary(name_file=fluidum_dir*"src/kernels/charm_reso_dubla_Kj/particles.data",output_path=fluidum_dir*"src/kernels/charm_reso_dubla_Kj/",format="jld2",kernel_type="total",Maxmass=10.1)
#keys(dic_charm)
export dic_charm
#keys(dic_NP)
#dic_NP.Jpsi.degeneracy