

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

function get_kernels(Tmin,Tmax,dT;type = "thermal", file_path=pwd()*"/src/FastReso_OC_kernels/", MC = false)
    particle_name=readdlm(file_path*"particles.data",comment_char='#',comments=true)[:,1]
    particle_MC=readdlm(file_path*"particles.data",comment_char='#',comments=true)[:,13]
    @show particle_MC

    Trange = collect(Int64,Tmin:dT:Tmax)
    kernel_dict = Dict()
    for (name, MC_code) in zip(particle_name, particle_MC)
        if MC == true
            kernel_path=[(file_path*"Kernels/PDGid_"*string(MC_code)*"_"*type*"_T0."*string(T)*"_Kj.out") for T in Trange]
        else
            kernel_path=[(file_path*"Kernels/"*name*"_"*type*"_T0."*string(T)*"_Kj.out") for T in Trange]
        end            
        kernel_dict[Symbol(name)] = kernel_path
    end
    
    return particle_name, kernel_dict
end


function save_kernels(;Tmin=1400,Tmax=1420,dT=1, type = "thermal", file_path=pwd()*"/src/FastReso_OC_kernels/", MC = false)    
    T = collect(Int64,Tmin:dT:Tmax)./10000 # in GeV
    particle_name, kernel_path = get_kernels(Tmin, Tmax, dT; type = type, file_path=file_path, MC = MC)
    
    for name in particle_name
        a = [readdlm(kernel_path[Symbol(name)][i], Float64,skipstart=1) for i in eachindex(T)]#same for all
        pt_list=reshape(a[1][:,1],(length(unique(a[1][:,1])),length(unique(a[1][:,2]))))
        ur_list=reshape(a[1][:,2],(length(unique(a[1][:,1])),length(unique(a[1][:,2]))))
        ur=ur_list[1,:]
        pt=pt_list[:,1]
        @cast a_new[i,j,k] := a[i][j,k] #recasting a
        itp = [interpolate((T,pt,ur),reshape(a_new[:,:,n],(length(T),length(pt),length(ur))),Gridded(Linear())) for n in 3:size(a[1])[2]]
        ext = extrapolate.(itp, Ref(Flat()))    
        #@show ext
        interp_kernel=kernels(ext...)
        @save file_path*name*"_"*type*".jld" interp_kernel
    end    
end



function save_kernels_diff(;Tmin=1400,Tmax=1420,dT=1, type = "thermal", file_path=pwd()*"/src/FastReso_OC_kernels/")    
    T = collect(Int64,Tmin:dT:Tmax)./10000 # in GeV
    particle_name, kernel_path = get_kernels(Tmin, Tmax, dT; type = type, file_path=file_path)
    
    for name in particle_name
        a = [readdlm(kernel_path[Symbol(name)][i], Float64,skipstart=1) for i in eachindex(T)]
        
        pt_list=reshape(a[1][:,1],(length(unique(a[1][:,1])),length(unique(a[1][:,2]))))
        ur_list=reshape(a[1][:,2],(length(unique(a[1][:,1])),length(unique(a[1][:,2]))))
        ur=ur_list[1,:]
        pt=pt_list[:,1]
        @cast a_new[i,j,k] := a[i][j,k] #recasting a
        itp = [interpolate((T,pt,ur),reshape(a_new[:,:,n],(length(T),length(pt),length(ur))),Gridded(Linear())) for n in 3:size(a[1])[2]]
        ext = extrapolate.(itp, Ref(Flat()))
        interp_kernel=kernels_diff(ext...)
        @save file_path*name*"_"*type*".jld" interp_kernel
    end    
end
   

function save_all_particles(;Tmin=1400,Tmax=1420,dT=10, type = "thermal")
    kernel_path = kernel_folder #personal path to the kernels folder
    
    save_kernels(;Tmin=Tmin,Tmax=Tmax,dT=dT,type = type, file_path=string(kernel_path,"/FastReso_kernels/"), MC = true)
    save_kernels_diff(;Tmin=Tmin,Tmax=Tmax,dT=dT,type = type, file_path=string(kernel_path,"/FastReso_OC_kernels/"))
    save_kernels_diff(;Tmin=Tmin,Tmax=Tmax,dT=dT,type = type, file_path=string(kernel_path,"/FastReso_HC_kernels/"))
end 

#uncomment only to save new kernels
#save_all_particles(Tmin=1400,Tmax=1420,dT=10,type = "thermal")
#save_all_particles(Tmin=1400,Tmax=1420,dT=10,type = "total")

