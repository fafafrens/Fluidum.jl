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
function initial_conditions(cent1,cent2;x,y,list=charm_eos.particle_list,gridpoints=500,rmax=30,rdrop=8,norm=3,diff=false)
    if diff==false
        return initialize_fields(x, cent1, cent2, list; tau0 = 0.4, gridpoints=500,rmax=30,norm_temp=1, norm_coll=1, exp_tail = true)
    elseif diff==true
        return initialize_fields(x, y, cent1, cent2, list; tau0 = 0.4, gridpoints=500,rmax=30,norm_temp=1, norm_coll=1, exp_tail = true)
    end
    
end

#HadronResonaceGas_ccbar(HadronResonaceGas().particle_list,3)
