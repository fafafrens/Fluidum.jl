using JLD2

using DelimitedFiles: readdlm

include("generating_kernels.jl")

struct particle_attribute{S,T,U,V}
    name::S
    mass::T
    degeneracy::U
    thermal_kernel_ext::V
    total_kernel_ext::V
end

const beauty_kernel_path = "C:\\Users\\feder\\.julia\\dev\\Fluidum\\src\\kernels\\beauty_reso_dubla_Kj\\"
const light_kernel_path = "C:\\Users\\feder\\.julia\\dev\\Fluidum\\src\\kernels\\FastReso_kernels\\"

function get_LF_Attributes(;file_path=light_kernel_path)
    data=readdlm(file_path*"light.txt",comment_char='#',comments=true)
    names    =convert.(String,data[:,1])
    mass     = convert.(Float64,data[:,2])
    spin     =convert.(Float64,data[:,4])

    thermal_kernel_path=[]
    total_kernel_path=[]
    
    protons_thermal=load(light_kernel_path*"proton_thermal.jld")["protons_thermal"]
    protons_total=load(light_kernel_path*"proton_total.jld"  )["protons_total"]
    kaons_thermal=load(light_kernel_path*"kaon_thermal.jld" )["kaons_thermal"]
    kaons_total=load(light_kernel_path*"kaon_total.jld" )["kaons_total"]
    pions_thermal=load(light_kernel_path*"pion_thermal.jld" )["pions_thermal"]
    pions_total=load(light_kernel_path*"pion_total.jld")["pions_total"]

    push!(thermal_kernel_path,protons_thermal,kaons_thermal,pions_thermal)
    push!(total_kernel_path,protons_total,kaons_total,pions_total)  
    
    
    return StructArray(particle_attribute.(
        names,
        mass,
        2*spin.+1,
        thermal_kernel_path,
        total_kernel_path
        ))
end


function get_diff_Attributes(;file_path=beauty_kernel_path,Maxmass=0.939,Minmass=0.13,condition=x->true)
    data=readdlm(file_path*"OpenBeautyParticleList_corrJS.txt",comment_char='#',comments=true)
    names    =convert.(String,data[:,1])
    mass     = convert.(Float64,data[:,2])
    spin     =convert.(Float64,data[:,4])

    thermal_kernel_path=[]
    total_kernel_path=[]
    

    B5279plu_thermal=load(beauty_kernel_path*"B5279plu_thermal.jld")["B5279plu_thermal"]
    B5324plu_thermal=load(beauty_kernel_path*"B5324plu_thermal.jld")["B5324plu_thermal"]
    Bs5366zer_thermal=load(beauty_kernel_path*"Bs5366zer_thermal.jld")["Bs5366zer_thermal"]
    Ups9460zer_thermal=load(beauty_kernel_path*"Ups9460zer_thermal.jld")["Ups9460zer_thermal"]
    Ups10000zer_thermal=load(beauty_kernel_path*"Ups10000zer_thermal.jld")["Ups10000zer_thermal"]
    B5279plu_total=load(beauty_kernel_path*"B5279plu_total.jld")["B5279plu_total"]
    B5324plu_total=load(beauty_kernel_path*"B5324plu_total.jld")["B5324plu_total"]
    Bs5366zer_total=load(beauty_kernel_path*"Bs5366zer_total.jld")["Bs5366zer_total"]
    Ups9460zer_total=load(beauty_kernel_path*"Ups9460zer_total.jld")["Ups9460zer_total"]
    Ups10000zer_total=load(beauty_kernel_path*"Ups10000zer_total.jld")["Ups10000zer_total"]

    push!(thermal_kernel_path,B5279plu_thermal,
    B5324plu_thermal,
    Bs5366zer_thermal,
    Ups9460zer_thermal,
    Ups10000zer_thermal)
    push!(total_kernel_path,B5279plu_total,
    B5324plu_total,
    Bs5366zer_total,
    Ups9460zer_total,
    Ups10000zer_total)  

    return StructArray(particle_attribute.(
        names,
        mass,
        2*spin.+1,
        thermal_kernel_path,
        total_kernel_path
        ))
end
    
function dictionary()
    LF = get_LF_Attributes()
    OB = get_diff_Attributes()
    Key_Tuple_LF=(:proton,:kaon,:pion)
    Key_Tuple_OB=(:B5279plu,
    :B5324plu,
    :Bs5366zer,
    :Ups9460zer,
    :Ups10000zer)#,:Lcplus#,:Xic,:Omc)
    
    Tuple1 = (; zip(Key_Tuple_LF, LF)...) 
    Tuple2 = (; zip(Key_Tuple_OB, OB)...) 
    #return Tuple1
    return merge(Tuple1, Tuple2)
end 

dic = dictionary();
export dic
