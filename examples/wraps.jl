using Fluidum

"""
    fluidproperties(eos;ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0,M=1.5)
Set up the fluid properties: equation of states and transport coefficients with relaxation time
# Arguments
- eos: choose among IdealQCD, FluiduMEoS, HadronResonaceGas, waleckacondition, LatticeQCD, HadronResonaceGasNew, HadronResonaceGas_ccbar
- output is a structure with fields eos, shear, bulk, diffusion 
# Options
- ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.2,M=1.5,diff=false
- if diff==true, includes HQ diffusion 
"""
function fluidproperties(eos;ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.2,M=1.5,diff=false)
    if diff==false
        Fluidum.FluidProperties(eos,QGPViscosity(ηs,Cs),SimpleBulkViscosity(ζs,Cζ),ZeroDiffusion())
    elseif diff==true
        Fluidum.FluidProperties(eos,QGPViscosity(ηs,Cs),SimpleBulkViscosity(ζs,Cζ),HQdiffusion(DsT,M))
    end
end

#dump(fluidproperties(FluiduMEoS();diff=true))
charm_eos = HadronResonaceGas()
"""
function initial_conditions(eos;gridpoints=500,rmax=30,rdrop=8,norm=3,diff=false)
returns DiscreteFields(disc,disc_fields,phi), containing the radial discretization, the fluid fields structure and the initialized fields
- if diff==true, includes HQ diffusion 
"""
function initial_conditions(cent1,cent2;x=artifact"ic_fluid",y=artifact"ncoll",list=charm_eos.particle_list,gridpoints=500,rmax=30,rdrop=8,norm=3,diff=false)
    if diff==false
        initialize_fields(x, cent1, cent2, list; tau0 = 0.4, gridpoints=500,rmax=30,norm_temp=1, norm_coll=1, exp_tail = true)
    elseif diff==true
        initialize_fields(x, y, cent1, cent2, list; tau0 = 0.4, gridpoints=500,rmax=30,norm_temp=1, norm_coll=1, exp_tail = true)
    end
    
end

#HadronResonaceGas_ccbar(HadronResonaceGas().particle_list,3)

#oneshoot(oned_visc_hydro_discrete,matrxi1d_visc_HQ!,fluidpropery1,phi1,tspan)

#check maybe the x and y can be called from the ic cent1 cent2
function runFluidum_heavy(eos,cent1,cent2,x::TabulatedTrento{A,B},y::TabulatedTrento{A,B};
    ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.2,M=1.5,
    tau0=0.4,Tfo=0.156,maxtime=30,
    gridpoints=500,rmax=30,norm_temp=1, norm_coll=1, exp_tail = true
    #particle_list=[:D0, :Dplus, :Dstarplus, :Dsplus, :jpsi,:Lcplus, :Xic, :Omc ]
    )
    fluidproperties_nodiff=Fluidum.FluidProperties(eos,QGPViscosity(ηs,Cs),SimpleBulkViscosity(ζs,Cζ),ZeroDiffusion())
    #    fields = initialize_fields(x, cent1, cent2, list; tau0 = tau0, gridpoints=gridpoints,rmax=rmax,norm_temp=norm_temp, norm_coll=norm_coll, exp_tail = exp_tail)
    fluidproperties_diff=Fluidum.FluidProperties(eos,QGPViscosity(ηs,Cs),SimpleBulkViscosity(ζs,Cζ),HQdiffusion(DsT,M))
    fields = initialize_fields(x, y, cent1, cent2, list; tau0 = tau0, gridpoints=gridpoints,rmax=rmax,norm_temp=norm_temp, norm_coll=norm_coll, exp_tail = exp_tail)
    
    disc,disc_fields,phi = @unpack fields
    tspan = (tau0,maxtime)
    #oneshoot(disc_fields,matrxi1d_visc_HQ!,fluidproperties,phi,tspan)
    fo_surface_diff=freeze_out_routine(disc_fields,matrxi1d_visc_HQ!,fluidproperties_diff,phi,tspan,Tfo=Tfo)
    fo_surface_nodiff=freeze_out_routine(disc_fields,matrxi1d_visc_HQ!,fluidproperties_nodiff,phi,tspan,Tfo=Tfo)
    return(fo_surface_nodiff,fo_surface_diff)
end

function runFluidum_light(eos,cent1,cent2,x::TabulatedTrento{A,B};
    ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.2,M=1.5,
    tau0=0.4,Tfo=0.156,maxtime=30,
    gridpoints=500,rmax=30,norm_temp=1, norm_coll=1, exp_tail = true
    #particle_list=[:D0, :Dplus, :Dstarplus, :Dsplus, :jpsi,:Lcplus, :Xic, :Omc ]
    )
    fluidproperties_nodiff=Fluidum.FluidProperties(eos,QGPViscosity(ηs,Cs),SimpleBulkViscosity(ζs,Cζ),ZeroDiffusion())
    #    fields = initialize_fields(x, cent1, cent2, list; tau0 = tau0, gridpoints=gridpoints,rmax=rmax,norm_temp=norm_temp, norm_coll=norm_coll, exp_tail = exp_tail)
    fields = initialize_fields(x, cent1, cent2, list; tau0 = tau0, gridpoints=gridpoints,rmax=rmax,norm_temp=norm_temp, norm_coll=norm_coll, exp_tail = exp_tail)
    
    disc,disc_fields,phi = @unpack fields
    tspan = (tau0,maxtime)
    #oneshoot(disc_fields,matrxi1d_visc_HQ!,fluidproperties,phi,tspan)
    fo_surface_nodiff=freeze_out_routine(disc_fields,matrxi1d_visc!,fluidproperties_nodiff,phi,tspan,Tfo=Tfo)
    return(fo_surface_nodiff)
end

function compute_observables_heavy(eos,cent1,cent2,x::TabulatedTrento{A,B},y::TabulatedTrento{A,B};
    ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.2,M=1.5,
    tau0=0.4,Tfo=0.156,maxtime=30,
    gridpoints=500,rmax=30,norm_temp=1, norm_coll=1, exp_tail = true,
    particle_list=[:D0, :Dplus, :Dstarplus, :Dsplus, :jpsi,:Lcplus, :Xic, :Omc ],
    pt_max=6,pt_min=0.)

    fo_surface_nodiff,fo_surface_diff=runFluidum_heavy(eos,cent1,cent2,x::TabulatedTrento{A,B},y::TabulatedTrento{A,B},
    ηs=ηs,Cs=Cs,ζs=ζs,Cζ=Cζ,DsT=DsT,M=M,
    tau0=tau0,Tfo=Tfo,maxtime=maxtime,
    gridpoints=gridpoints,rmax=rmax,norm_temp=norm_temp, norm_coll=norm_coll, exp_tail = exp_tail
    )

    map(particle_list) do part
         
        g_th= spectra(fo_surface_nodiff,part,pt_max=pt_max,pt_min=pt_min,decays=false)
        g1_th= spectra(fo_surface_diff,part,pt_max=pt_max,pt_min=pt_min,decays=false)

        g_tot= spectra(fo_surface_nodiff,part,pt_max=pt_max,pt_min=pt_min,decays=true)
        g1_tot= spectra(fo_surface_diff,part,pt_max=pt_max,pt_min=pt_min,decays=true)
        
        x = range(pt_min,pt_max,100)

        open("dNdptpt"*string(part)*"_"*string(cent1)*"_"*string(cent2)*"_"*string(Tfo)*".txt", "w") do io
            write(io, "# pt\t th_nodiff\t th_diff\t dec_nodiff\t dec_diff\n")
            writedlm(io, [x g_th g1_th g_tot g1_tot])
        end
    end

    open("dNdy_"*"_"*string(cent1)*"_"*string(cent2)*"_"*string(Tfo)*".txt", "w") do io
        write(io, "# part_name \t th_int_yield \t dec_int_yield\n")
           
        map(particle_list) do part
            yield_th, err=multiplicity(fo_surface_nodiff,dic[part];decays=false)
            yield_tot, err=multiplicity(fo_surface_nodiff,dic[part];decays=true)
            writedlm(io, [string(part)  yield_th  yield_tot])            
        end
    end    


end
  