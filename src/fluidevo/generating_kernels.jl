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



function get_kernels(T;type = "thermal", file_path=pwd()*"/src/FastReso_OC_kernels/", MC = false)
    particle_name=readdlm(file_path*"particles.data",comment_char='#',comments=true)[:,1]
    particle_MC=readdlm(file_path*"particles.data",comment_char='#',comments=true)[:,13]
    
    
    kernel_dict = Dict()
    for (name, MC_code) in zip(particle_name, particle_MC)
        if MC == true
            kernel_path=(file_path*"Kernels/PDGid_"*string(MC_code)*"_"*type*"_T0."*string(T)*"_Kj.out") 
        else
            kernel_path=(file_path*"Kernels/"*name*"_"*type*"_T0."*string(T)*"_Kj.out") 
        end            
        kernel_dict[Symbol(name)] = kernel_path
    end
    
    return particle_name, kernel_dict
end


function save_kernels(;T=1100,type = "thermal", file_path=pwd()*"/src/FastReso_OC_kernels/", MC = false)    
    particle_name, kernel_path = get_kernels(T; type = type, file_path=file_path, MC = MC)
    
    for name in particle_name
        a = readdlm(kernel_path[Symbol(name)], Float64,skipstart=1) #same for all
        pt_list=reshape(a[:,1],(length(unique(a[:,1])),length(unique(a[:,2])))) #riformatting so that pt_list is a matrix with all equal columns and rows with pt entries
        ur_list=reshape(a[:,2],(length(unique(a[:,1])),length(unique(a[:,2])))) #riformatting so that ur_list is a matrix with all equal rows and columns with ur entries
        ur=ur_list[1,:]
        pt=pt_list[:,1]
        
        itp = [interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear())) for n in 3:size(a)[2]] #interpolation of all the columns from 3 on, through linear interpolation
        ext = extrapolate.(itp, Ref(Flat()))    
        
        interp_kernel=kernels(ext...) #performing "unpaking": taking all the elements of ext, and using them to fill kernels 
        @save file_path*name*"_"*type*".jld" interp_kernel
    end    
end

function save_kernels_diff(;T=1100,type = "thermal", file_path=pwd()*"/src/FastReso_OC_kernels/", MC = false)    
    particle_name, kernel_path = get_kernels(T; type = type, file_path=file_path, MC = MC)
    
    for name in particle_name
        a = readdlm(kernel_path[Symbol(name)], Float64,skipstart=1) #same for all
        pt_list=reshape(a[:,1],(length(unique(a[:,1])),length(unique(a[:,2])))) 
        ur_list=reshape(a[:,2],(length(unique(a[:,1])),length(unique(a[:,2])))) 
        ur=ur_list[1,:]
        pt=pt_list[:,1]
        
        itp = [interpolate((pt,ur),reshape(a[:,n],(length(pt),length(ur))),Gridded(Linear())) for n in 3:size(a)[2]] 
        ext = extrapolate.(itp, Ref(Flat()))    
        interp_kernel=kernels_diff(ext...) 
        @save file_path*name*"_"*type*".jld" interp_kernel
    end    
end


function save_all_particles(;T=1100, type = "thermal")
    kernel_path = kernel_folder #personal path to the kernels folder
    
    save_kernels(;T=T,type = type, file_path=string(kernel_path,"/FastReso_kernels/"), MC = true)
    save_kernels_diff(;T=T,type = type, file_path=string(kernel_path,"/FastReso_OC_kernels/"), MC = false)
    save_kernels_diff(;T=T,type = type, file_path=string(kernel_path,"/FastReso_HC_kernels/"), MC = false)
end 


#uncomment only to save new kernels
#save_all_particles(;T=1500, type = "thermal")
#save_all_particles(;T=1500, type = "total")


