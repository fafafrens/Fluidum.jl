


struct particle_attribute{S,T,U,V}
    name::S
    mass::T
    degeneracy::U
    thermal_kernel_ext::V
    total_kernel_ext::V
end



function load_kernel(file_path,name;type="thermal")
    return load(file_path*name*"_"*type*".jld"; typemap=Dict("Main.kernels_diff" => kernels_diff))["interp_kernel"]
    #return load(file_path*name*"_"*type*".jld"; typemap=Dict("Main.kernels_diff" => kernels_diff))["interp_kernel"]
end


function get_LF_Attributes(;file_path=pwd()*"/src/FastReso_kernels/")
    data=readdlm(file_path*"particles.data",comment_char='#',comments=true)
    names    =convert.(String,data[:,1])
    mass     = convert.(Float64,data[:,2])
    spin     =convert.(Float64,data[:,4])

    thermal_kernel_path=[]
    total_kernel_path=[]
    
    # @load "C:\\Data\\Università\\Germania\\Tesi_mag\\solver\\src\\FieldsEvolution\\src\\Kernels_Folder/FastReso_kernels/proton_thermal.jld"
    # @load "C:\\Data\\Università\\Germania\\Tesi_mag\\solver\\src\\FieldsEvolution\\src\\Kernels_Folder/FastReso_kernels/kaon_thermal.jld"
    # @load "C:\\Data\\Università\\Germania\\Tesi_mag\\solver\\src\\FieldsEvolution\\src\\Kernels_Folder/FastReso_kernels/pion_thermal.jld"

    # @load "C:\\Data\\Università\\Germania\\Tesi_mag\\solver\\src\\FieldsEvolution\\src\\Kernels_Folder/FastReso_kernels/proton_total.jld"
    # @load "C:\\Data\\Università\\Germania\\Tesi_mag\\solver\\src\\FieldsEvolution\\src\\Kernels_Folder/FastReso_kernels/kaon_total.jld"
    # @load "C:\\Data\\Università\\Germania\\Tesi_mag\\solver\\src\\FieldsEvolution\\src\\Kernels_Folder/FastReso_kernels/pion_total.jld" 

    # push!(thermal_kernel_path,protons_thermal,kaons_thermal,pions_thermal)
    
    # push!(total_kernel_path,protons_total,kaons_total,pions_total)  
    for name in names 
        thermal = load(file_path*name*"_thermal.jld"; typemap=Dict("Main.kernels" => kernels))["interp_kernel"]
        total = load(file_path*name*"_total.jld"; typemap=Dict("Main.kernels" => kernels))["interp_kernel"]
        push!(thermal_kernel_path,thermal)
        push!(total_kernel_path,total)  
    end
    
    return StructArray(particle_attribute.(
        names,
        mass,
        2*spin.+1,
        thermal_kernel_path,
        total_kernel_path
        ))
end


function get_diff_Attributes(;file_path=pwd()*"/src/FastReso_OC_kernels/",Maxmass=0.939,Minmass=0.13,condition=x->true)
    data=readdlm(file_path*"particles.data",comment_char='#',comments=true)
    names    =convert.(String,data[:,1])
    mass     = convert.(Float64,data[:,2])
    spin     =convert.(Float64,data[:,4])

    thermal_kernel_path=[]
    total_kernel_path=[]
    for name in names 
        thermal = load(file_path*name*"_thermal.jld"; typemap=Dict("Main.kernels_diff" => kernels_diff))["interp_kernel"]
        total = load(file_path*name*"_total.jld"; typemap=Dict("Main.kernels_diff" => kernels_diff))["interp_kernel"]
        push!(thermal_kernel_path,thermal)
        push!(total_kernel_path,total)  
    end
    
    return StructArray(particle_attribute.(
        names,
        mass,
        2*spin.+1,
        thermal_kernel_path,
        total_kernel_path
        ))
end
    
function dictionary()
    LF = get_LF_Attributes(;file_path="$(@__DIR__)/Kernels_folder/FastReso_kernels/")
    OC = get_diff_Attributes(;file_path="$(@__DIR__)/Kernels_folder/FastReso_OC_kernels/")
    HC = get_diff_Attributes(;file_path="$(@__DIR__)/Kernels_folder/FastReso_HC_kernels/")[1]
    Key_Tuple_LF=(:proton,:kaon,:pion)
    Key_Tuple_OC=(:D0,:Dplus,:Dstar0,:Dstarplus,:Dsplus)#,:Lcplus#,:Xic,:Omc)
    Key_Tuple_HC=(:Jpsi)
    
    Tuple1 = (; zip(Key_Tuple_LF, LF)...) 
    Tuple2 = (; zip(Key_Tuple_OC, OC)...) 
    Tuple3 = (jpsi = HC,) 
    #return Tuple1
    return merge(Tuple1, Tuple2, Tuple3)
end 

dic = dictionary();
export dic
