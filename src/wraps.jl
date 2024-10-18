abstract type ProfilePars end

struct StdPars <: ProfilePars
    σ::Float64
    norm::Float64
end

StdPars(;σ=1,norm=1) = StdPars(σ,norm)
struct TrentoPars <: ProfilePars
    cent1::Int64
    cent2::Int64
    norm::Float64
    filepath::String
end

TrentoPars(;norm=1,filepath=pwd()*"/../examples/ic_data/only_shift_BKG_full_order_changed.txt",cent1=0,cent2=10) = TrentoPars(cent1,cent2,norm,filepath)

struct InterpolPars <: ProfilePars
    norm::Float64
    filepath::String
end

InterpolPars(;norm=1,filepath=pwd()*"/../examples/ic_data/tempBG502_0-10.txt") = InterpolPars(norm,filepath)

struct FugacityPars <: ProfilePars
    norm_coll::Float64 
    σ_in::Float64
    dσ_QQdy::Float64
    rdrop::Float64  
    n0::Float64
end

FugacityPars(;norm_coll=1, σ_in=70,dσ_QQdy=0.463,rdrop=4, n0=31.57) = FugacityPars(norm_coll,σ_in, dσ_QQdy, rdrop, n0)

#chiamale tutte uguali ma con tipo di parametro di tipo diverso
function gauss(x;pars=StdPars())
    exp(-x^2/(2pars.σ^2))*pars.norm+0.1
end

function fermidirac(x;pars=StdPars())
    return 1/(1+exp(pars.σ*x))*pars.norm+0.1
end

function trento_profile(x;pars=TrentoPars())
    return Profiles(TabulatedTrento(pars.filepath),pars.cent1,pars.cent2,norm_temp=pars.norm)(x)
end

function temp_interpol(x;pars=InterpolPars())
    file=readdlm(pars.filepath)
    #TabulatedTrento(file[:,1], file[:,2:end])
    pars.norm.*linear_interpolation(reverse(file[:,1]),reverse(file[:,2]); extrapolation_bc=Flat())(x).+0.0001
end

function ncoll(r,pars)
    #rdrop=4,n0 = 31.57)
    #return 7*n0*(initial_temp(r)/initial_temp(0.0))^4/3 +0.001
    if r<=pars.rdrop
        return pars.n0
    else 
        return 0.0001 #fm-2 since n0 [fm-2]
    end
end

function nhard(x,pars::FugacityPars,temp_pars::ProfilePars,temp_profile::Symbol,tau0)
    temp(x)=profile_functions[temp_profile](x,pars=temp_pars)
    pars.norm_coll*2/tau0/pars.σ_in*pars.dσ_QQdy*ncoll(0,pars)*(temp(x)/temp(0.))^4 +0.001
end 

function fugacity(x;pars=FugacityPars(),temp_profile::Symbol=:t_int,temp_pars=InterpolPars(),tau0=0.4,eos=HadronResonaceGas())
    temp(x)=profile_functions[temp_profile](x,pars=temp_pars)
    if x<=pars.rdrop
        return log(nhard(x,pars,temp_pars,temp_profile,tau0)/(thermodynamic(temp(x),0.0,eos).pressure)).+ 0.0001
        else return log(nhard(pars.rdrop,pars,temp_pars,temp_profile,tau0)/(thermodynamic(temp(pars.rdrop),0.0,eos).pressure)).+ 0.0001
        end
end 
# Dictionary to map profile names to functions
const profile_functions = Dict(
    :gauss => gauss,
    :fd => fermidirac,
    :trento => trento_profile,
    :t_int => temp_interpol,
    :fug => fugacity
)

function check(profile)
    if !haskey(profile_functions, profile)
        error("Unknown profile: $profile")
    end
end

function apply_profile(profile_name::Symbol, x; pars=nothing)
    if !haskey(profile_functions, profile_name)
        error("Unknown profile: $profile_name")
    end

    profile_func = profile_functions[profile_name]
   
    # Apply the profile function with parameters, using default if params is nothing
    if params === nothing
        return profile_func(x)
    else
        return profile_func(x; pars=pars)
    end
end

function initial_conditions(;eos_HQ=HadronResonaceGas(),gridpoints=500,rmax=30,
    temp_profile::Symbol=:t_int, fug_profile::Symbol=:fug, fug_pars=FugacityPars(), temp_pars=InterpolPars(),
    tau0=0.4)


    check(temp_profile)
    check(fug_profile)
    
    temp(x)=profile_functions[temp_profile](x,pars=temp_pars)
    fug(x) = profile_functions[fug_profile](x,pars=fug_pars,temp_profile=temp_profile,temp_pars=temp_pars,tau0=tau0,eos=eos_HQ)
    
    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    oned_visc_hydro = Fluidum.HQ_viscous_1d()
    disc_fields = DiscreteFileds(oned_visc_hydro,disc,Float64) 
    phi=set_array((x)->temp(x),:temperature,disc_fields); #temperature initialization
    set_array!(phi,(x)->fug(x),:mu,disc_fields); #fugacity initialization
    NQQ̄,err= quadgk(x->2*pi*x*tau0*thermodynamic(temp(x),fug(x),eos_HQ).pressure,0,rmax,rtol=0.00001)
    @show NQQ̄
    
    #@show phi
return DiscreteFields(disc,disc_fields,phi)
end

#**************ic data la centralità
#**************ic dato il percorso file


#give kwargs for a generic profile function as additional optional arguments 

function runFluidum_fo(eos;eos_HQ=HadronResonaceGas(),
    ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.2,M=1.5,
    tau0=0.4,Tfo=0.156,maxtime=30,
    gridpoints=500,rmax=30,
    temp_profile::Symbol=:t_int, fug_profile::Symbol=:fug, fug_pars=FugacityPars(), temp_pars=InterpolPars())
    if iszero(DsT)
        fluidproperties=Fluidum.FluidProperties(eos,QGPViscosity(ηs,Cs),SimpleBulkViscosity(ζs,Cζ),ZeroDiffusion())
    else
        fluidproperties=Fluidum.FluidProperties(eos,QGPViscosity(ηs,Cs),SimpleBulkViscosity(ζs,Cζ),HQdiffusion(DsT,M))
    end

    #@show thermodynamic(temp_profile(0,norm_temp,σ=σ_temp),0.0,eos_HQ).pressure
    fields=initial_conditions(;eos_HQ=eos_HQ,gridpoints=gridpoints,rmax=rmax,
    temp_profile=temp_profile, fug_profile=fug_profile, fug_pars=fug_pars, temp_pars=temp_pars,
    tau0=tau0)
    tspan = (tau0,maxtime)
    if fields.initial_field[1,1]<Tfo
        println("Error: Tfo = "*string(Tfo)*" MeV is larger than max temperature in the inital profile T0 = "*string(fields.initial_field[1,1]))
        return nothing
    end
    return (fo=freeze_out_routine(fields.discrete_field,Fluidum.matrxi1d_visc_HQ!,fluidproperties,fields.initial_field,tspan,Tfo=Tfo),fluidproperties=fluidproperties)
end

function runFluidum(eos;eos_HQ=HadronResonaceGas(),
    ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.2,M=1.5,
    tau0=0.4,maxtime=30,
    gridpoints=500,rmax=30,
    temp_profile::Symbol=:t_int, fug_profile::Symbol=:fug, fug_pars=FugacityPars(), temp_pars=InterpolPars())
    if DsT == 0
        fluidproperties=Fluidum.FluidProperties(eos,QGPViscosity(ηs,Cs),SimpleBulkViscosity(ζs,Cζ),ZeroDiffusion())
       # norm_coll=0.
    else
        fluidproperties=Fluidum.FluidProperties(eos,QGPViscosity(ηs,Cs),SimpleBulkViscosity(ζs,Cζ),HQdiffusion(DsT,M))
    end
    fields=initial_conditions(;eos_HQ=eos_HQ,gridpoints=gridpoints,rmax=rmax,
    temp_profile=temp_profile, fug_profile=fug_profile, fug_pars=fug_pars, temp_pars=temp_pars,
    tau0=tau0)
    tspan = (tau0,maxtime)
    return oneshoot(fields.discrete_field,Fluidum.matrxi1d_visc_HQ!,fluidproperties,fields.initial_field,tspan)

end

function fields_evolution(eos;eos_HQ=HadronResonaceGas(),
    ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.2,M=1.5,
    tau0=0.4,Tfo=0.156,maxtime=30,
    gridpoints=500,rmax=30,
    temp_profile::Symbol=:t_int, fug_profile::Symbol=:fug, fug_pars=FugacityPars(), temp_pars=InterpolPars())
    
    fullevo=runFluidum(eos;eos_HQ=eos_HQ,
    ηs=ηs,Cs=Cs,ζs=ζs,Cζ=Cζ,DsT=DsT,M=M,
    tau0=tau0,maxtime=maxtime,
    gridpoints=gridpoints,rmax=rmax,
    temp_profile=temp_profile, fug_profile=fug_profile, fug_pars= fug_pars, temp_pars=temp_pars)
    #dump(Fields(fullevo))
    return Fields(fullevo)
    
end

struct Observables{S,T,U,V,M,K,A,B,C,D}
    particle::particle_attribute{S,T,U,V}
    yield_th::Float64
    yield_tot::Float64
    pt_bins::M
    spectra_th::K
    spectra_tot::K
    fluid_properties::FluidProperties{A,B,C,D}
    Tfo::Float64
end

function Observables(fo::FreezeOutResult{M,N},part::particle_attribute{S,T,U,V},fluidproperties::FluidProperties{A,B,C,D},Tfo;
    pt_min=0,pt_max=10.0,step=100,ccbar=0.7) where {M,N,S,T,U,V,A,B,C,D}
    spectra_th=getindex.(spectra(fo,part,pt_max=pt_max,pt_min=pt_min,step=step,decays=false,ccbar=ccbar)[:],1)
    #@show spectra(fo,part,pt_max=pt_max,pt_min=pt_min,step=step,decays=false,ccbar=ccbar)
    spectra_tot=getindex.(spectra(fo,part,pt_max=pt_max,pt_min=pt_min,step=step,decays=true,ccbar=ccbar)[:],1)
    yield_th, err=multiplicity(fo,part,decays=false)
    yield_tot, err=multiplicity(fo,part,decays=true)
    pt_bins = range(pt_min,pt_max,step) 
    return Observables(part,yield_th,yield_tot,pt_bins,spectra_th,spectra_tot,fluidproperties,Tfo)
end

function Observables(fo::FreezeOutResult{M,N},m::Float64,fluidproperties::FluidProperties{A,B,C,D},Tfo;
    pt_min=0,pt_max=10.0,step=100,deg=1, charge = 1) where {M,N,A,B,C,D}
    spectra_th=getindex.(spectra_internal(m,fo,pt_max=pt_max,pt_min=pt_min,step=step,deg=deg)[:],1)
    spectra_tot=spectra_th
    yield_th, err=multiplicity(m,fo)
    yield_tot = yield_th
    pt_bins = range(pt_min,pt_max,step) 
    part = particle_attribute(string(m),m,deg,charge, nothing,nothing)
    return Observables(part,yield_th,yield_tot,pt_bins,spectra_th,spectra_tot,fluidproperties,Tfo)
end

function compute_observables(eos,part::particle_attribute{S,T,U,V};eos_HQ=HadronResonaceGas(),
    ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.2,M=1.5,
    tau0=0.4,Tfo=0.156,maxtime=30,
    gridpoints=500,rmax=30,
    ccbar = 1.4,
    temp_profile::Symbol=:t_int, fug_profile::Symbol=:fug, fug_pars=FugacityPars(), temp_pars=InterpolPars(),
    pt_min=0,pt_max=10.0,step=100,save=false,savepath="./results/") where {S,T,U,V}

    fo, fluidproperties=runFluidum_fo(eos;eos_HQ=eos_HQ,
    ηs=ηs,Cs=Cs,ζs=ζs,Cζ=Cζ,DsT=DsT,M=M,
    tau0=tau0,Tfo=Tfo,maxtime=maxtime,
    gridpoints=gridpoints,rmax=rmax,
    temp_profile=temp_profile, fug_profile=fug_profile, fug_pars= fug_pars, temp_pars=temp_pars)

    obs = Observables(fo,part,fluidproperties, Tfo,pt_min=pt_min,pt_max=pt_max,step=step,ccbar = ccbar)

    if save==true
        save_observables(obs,path=savepath)
    end

    return obs
end

function compute_observables(eos,m::Float64;eos_HQ=HadronResonaceGas(),
    ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.2,M=1.5,
    tau0=0.4,Tfo=0.156,maxtime=30,
    gridpoints=500,rmax=30,
    temp_profile::Symbol=:t_int, fug_profile::Symbol=:fug, fug_pars=FugacityPars(), temp_pars=InterpolPars(),
    pt_min=0,pt_max=10.0,step=100,save=false,savepath="./results/") 

    fo, fluidproperties=runFluidum_fo(eos;eos_HQ=eos_HQ,
    ηs=ηs,Cs=Cs,ζs=ζs,Cζ=Cζ,DsT=DsT,M=M,
    tau0=tau0,Tfo=Tfo,maxtime=maxtime,
    gridpoints=gridpoints,rmax=rmax,
    temp_profile=temp_profile, fug_profile=fug_profile, fug_pars= fug_pars, temp_pars=temp_pars)

    obs = Observables(fo,m,fluidproperties,Tfo, pt_min=pt_min,pt_max=pt_max,step=step)

    if save==true
        save_observables(obs,path=savepath)
    end

    return obs
end


function save_observables(obj::Observables{S,T,U,V,M,K,A,B,C,D};path="./results/") where {S,T,U,V,M,K,A,B,C,D}
    

    if isdir(path)
    else mkdir(path)
    end

    filename =  get_filename(obj,path=path)
    @show filename
    #if isfile(filename)
    #    println("file already exists")
    #    return nothing
    #else
    open(filename, "w") do io 
        write(io, "# int_yield_th: "*string(obj.yield_th)*", int_yield_tot: "*string(obj.yield_tot)*"\n")
        write(io, "# pt\t th\t tot\t \n")
        writedlm(io, [obj.pt_bins obj.spectra_th obj.spectra_tot])
        close(io)
    #end
end
end

function get_filename(obj::Observables{S,T,U,V,M,K,A,B,C,D};path="./results/",toplot=false) where {S,T,U,V,M,K,A,B,C,D}
    
    part = obj.particle.name
    if typeof(obj.fluid_properties.shear)<:ZeroViscosity
        ηs = 0
    else    
        ηs = obj.fluid_properties.shear.ηs
    end
    if typeof(obj.fluid_properties.bulk)<:ZeroBulkViscosity
        ζs = 0
    else    
        ζs = obj.fluid_properties.bulk.ζs
    end
    if typeof(obj.fluid_properties.diffusion)<:ZeroDiffusion
        DsT = 0
    else    
        DsT = obj.fluid_properties.diffusion.DsT
    end
    Tfo = obj.Tfo

    if toplot == true
        return path*string(part)*"_Tfo_"*string(Tfo)*"_ηs_"*string(ηs)*"_ζs_"*string(ζs)*"_DsT_"*string(DsT)*"_observables.pdf"
    else
        return path*string(part)*"_Tfo_"*string(Tfo)*"_ηs_"*string(ηs)*"_ζs_"*string(ζs)*"_DsT_"*string(DsT)*"_observables.txt"
    end
end

function read_data(file)
    readdlm(file)
end

