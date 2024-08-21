function gauss(x,norm;σ=1)
    return exp(-x^2/(2σ^2))*norm+0.1
end

function fermidirac(x,norm;σ=1)
    return 1/(1+exp(σ*x))*norm+0.1
end

function trento_profile(x,norm;σ=1,filepath=pwd()*"/../examples/ic_data/only_shift_BKG_full_order_changed.txt",cent1=0,cent2=10)
    return Profiles(TabulatedTrento(filepath),cent1,cent2,norm_temp=norm)(x)
end

function initial_conditions(;gridpoints=500,rmax=30,
    norm_temp=0.5, σ_temp=1, norm_coll=-4, σ_fug=1,temp_profile=gauss, fug_profile=gauss)
    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    oned_visc_hydro = Fluidum.HQ_viscous_1d()
    disc_fields = DiscreteFileds(oned_visc_hydro,disc,Float64) 
    phi=set_array((x)->temp_profile(x,norm_temp;σ=σ_temp),:temperature,disc_fields); #temperature initialization
    set_array!(phi,(x)->fug_profile(x,norm_coll;σ=σ_fug),:mu,disc_fields); #fugacity initialization
return DiscreteFields(disc,disc_fields,phi)
end

#**************ic data la centralità
#**************ic dato il percorso file

function runFluidum_fo(eos;
    ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.2,M=1.5,
    tau0=0.4,Tfo=0.156,maxtime=30,
    gridpoints=500,rmax=30,
    norm_temp=0.5, σ_temp=1, norm_coll=-4, σ_fug=1,temp_profile=gauss,fug_profile=gauss)
    if DsT == 0
        fluidproperties=Fluidum.FluidProperties(eos,QGPViscosity(ηs,Cs),SimpleBulkViscosity(ζs,Cζ),ZeroDiffusion())
        norm_coll=0.
    else
        fluidproperties=Fluidum.FluidProperties(eos,QGPViscosity(ηs,Cs),SimpleBulkViscosity(ζs,Cζ),HQdiffusion(DsT,M))
    end
    fields=initial_conditions(;gridpoints=gridpoints,rmax=rmax,norm_temp=norm_temp, σ_temp=σ_temp, norm_coll=norm_coll, σ_fug=σ_fug,temp_profile=temp_profile,fug_profile=fug_profile)
    tspan = (tau0,maxtime)
    if fields.initial_field[1,1]<Tfo
        println("Error: Tfo = "*string(Tfo)*" MeV is larger than max temperature in the inital profile T0 = "*string(fields.initial_field[1,1]))
        return nothing
    end
    return (fo=freeze_out_routine(fields.discrete_field,Fluidum.matrxi1d_visc_HQ!,fluidproperties,fields.initial_field,tspan,Tfo=Tfo),fluidproperties=fluidproperties)
end


function runFluidum(eos;
    ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.2,M=1.5,
    tau0=0.4,maxtime=30,
    gridpoints=500,rmax=30,
    norm_temp=0.5, σ_temp=1, norm_coll=-4, σ_fug=1,temp_profile=gauss,fug_profile=gauss)
    if DsT == 0
        fluidproperties=Fluidum.FluidProperties(eos,QGPViscosity(ηs,Cs),SimpleBulkViscosity(ζs,Cζ),ZeroDiffusion())
        norm_coll=0.
    else
        fluidproperties=Fluidum.FluidProperties(eos,QGPViscosity(ηs,Cs),SimpleBulkViscosity(ζs,Cζ),HQdiffusion(DsT,M))
    end
    fields=initial_conditions(;gridpoints=gridpoints,rmax=rmax,norm_temp=norm_temp, σ_temp=σ_temp, norm_coll=norm_coll, σ_fug=σ_fug,temp_profile=temp_profile,fug_profile=fug_profile)
    tspan = (tau0,maxtime)
    return oneshoot(fields.discrete_field,Fluidum.matrxi1d_visc_HQ!,fluidproperties,fields.initial_field,tspan)

end

function fields_evolution(eos;
    ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.2,M=1.5,
    tau0=0.4,maxtime=30,
    gridpoints=500,rmax=30,
    norm_temp=0.5, σ_temp=1, norm_coll=-4, σ_fug=1,temp_profile=gauss,fug_profile=gauss)
    if DsT == 0
        norm_coll=0
    end
    fullevo=runFluidum(eos,ηs=ηs,Cs=Cs,ζs=ζs,Cζ=Cζ,DsT=DsT,M=M,
    tau0=tau0,maxtime=maxtime,
    gridpoints=gridpoints,rmax=rmax,
    norm_temp=norm_temp, σ_temp=σ_temp, norm_coll=norm_coll, σ_fug=σ_fug,temp_profile=temp_profile,fug_profile=fug_profile)
    #dump(Fields(fullevo))
    #return Fields(fullevo)
    
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
    pt_min=0,pt_max=10.0,step=100) where {M,N,S,T,U,V,A,B,C,D}
    spectra_th,err=spectra(fo,part,pt_max=pt_max,pt_min=pt_min,step=step,decays=false)
    spectra_tot,err=spectra(fo,part,pt_max=pt_max,pt_min=pt_min,step=step,decays=true)
    yield_th, err=multiplicity(fo,part,decays=false)
    yield_tot, err=multiplicity(fo,part,decays=true)
    pt_bins = range(pt_min,pt_max,step) 
    return Observables(part,yield_th,yield_tot,pt_bins,spectra_th,spectra_tot,fluidproperties,Tfo)
end

function Observables(fo::FreezeOutResult{M,N},m::Float64,fluidproperties::FluidProperties{A,B,C,D},Tfo;
    pt_min=0,pt_max=10.0,step=100,deg=1) where {M,N,A,B,C,D}
    spectra_th=getindex.(spectra_internal(m,fo,pt_max=pt_max,pt_min=pt_min,step=step,deg=deg)[:],1)
    spectra_tot=spectra_th
    yield_th, err=multiplicity(m,fo)
    yield_tot = yield_th
    pt_bins = range(pt_min,pt_max,step) 
    part = particle_attribute(string(m),m,deg,nothing,nothing)
    return Observables(part,yield_th,yield_tot,pt_bins,spectra_th,spectra_tot,fluidproperties,Tfo)
end

function compute_observables(eos,part::particle_attribute{S,T,U,V};
    ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.2,M=1.5,
    tau0=0.4,Tfo=0.156,maxtime=30,
    gridpoints=500,rmax=30,
    norm_temp=0.5, σ_temp=1, norm_coll=-4, σ_fug=1,temp_profile=gauss,fug_profile=gauss,
    pt_min=0,pt_max=10.0,step=100,save=false,savepath="./results/") where {S,T,U,V}

    fo, fluidproperties=runFluidum_fo(eos,ηs=ηs,Cs=Cs,ζs=ζs,Cζ=Cζ,DsT=DsT,M=M,
    tau0=tau0,Tfo=Tfo,maxtime=maxtime,
    gridpoints=gridpoints,rmax=rmax,
    norm_temp=norm_temp, σ_temp=σ_temp, norm_coll=norm_coll, σ_fug=σ_fug,temp_profile=temp_profile,fug_profile=fug_profile)

    obs = Observables(fo,part,fluidproperties, pt_min=pt_min,pt_max=pt_max,step=step,Tfo)

    if save==true
        save_observables(obs,path=savepath)
    end

    return obs
end

function compute_observables(eos,m::Float64;
    ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.2,M=1.5,
    tau0=0.4,Tfo=0.156,maxtime=30,
    gridpoints=500,rmax=30,
    norm_temp=0.5, σ_temp=1, norm_coll=-4, σ_fug=1,temp_profile=gauss,fug_profile=gauss,
    pt_min=0,pt_max=10.0,step=100,save=false,savepath="./results/") where {S,T,U,V}

    fo, fluidproperties=runFluidum_fo(eos,ηs=ηs,Cs=Cs,ζs=ζs,Cζ=Cζ,DsT=DsT,M=M,
    tau0=tau0,Tfo=Tfo,maxtime=maxtime,
    gridpoints=gridpoints,rmax=rmax,
    norm_temp=norm_temp, σ_temp=σ_temp, norm_coll=norm_coll, σ_fug=σ_fug,temp_profile=temp_profile,fug_profile=fug_profile)

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
    if isfile(filename)
        println("file already exists")
        return nothing
    else
    open(filename, "w") do io 
        write(io, "# int_yield_th: "*string(obj.yield_th)*", int_yield_tot: "*string(obj.yield_tot)*"\n")
        write(io, "# pt\t th\t tot\t \n")
        writedlm(io, [obj.pt_bins obj.spectra_th obj.spectra_tot])
        close(io)
    end
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

