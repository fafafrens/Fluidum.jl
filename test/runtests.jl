using Fluidum
using Test
using LinearAlgebra
using ForwardDiff
using QuadGK
using DifferentiationInterface


@testset "equation of state" begin
    ccbar = 30.
    @test IdealQCD()==IdealQCD(3,2)
    @test isfinite(thermodynamic(1,IdealQCD()).pressure[1])
    @test isfinite(thermodynamic(0,IdealQCD()).pressure_derivative[1])
  
    @test length(thermodynamic(10,IdealQCD(1,0)).pressure)==1
    @test length(thermodynamic(10,IdealQCD(1,0)).pressure_derivative)==1
    @test length(thermodynamic(10,IdealQCD(1,0)).pressure_hessian)==1

    @test isfinite(free_charm(1,0,Heavy_Quark(readresonancelist(), ccbar))[1])
    @test round(pressure(1,Heavy_Quark(readresonancelist(), ccbar)), sigdigits=10)==round(pressure(1,FluiduMEoS()), sigdigits=10)
    @test round(pressure(1,Heavy_Quark(readresonancelist(), ccbar)), sigdigits=10)==round(thermodynamic(1,FluiduMEoS()).pressure[1], sigdigits=10)
    @test pressure(1,Heavy_Quark(readresonancelist(), ccbar)) ≈ pressure(1,FluiduMEoS()) atol=0.01

    eos=FluiduMEoS()
    @test  isapprox(DifferentiationInterface.derivative(x->Fluidum.pressure_derivative(x,Val(2),eos),AutoForwardDiff(),0.2),Fluidum.pressure_derivative(0.2,Val(3),eos))
end


@testset "fluid_properties" begin
    eos = FluiduMEoS()
    @test det(Fluidum.one_d_viscous_matrix_fugacity([0.2,0.1,0,0,0,-0.1,0],2.,2.,0.5,0.5,0,0,0,0,0,0,1.,1.,0.1,0.1,0.1,1.,0)[1])!=0 
    params=Fluidum.FluidProperties(eos,QGPViscosity(0.2,0.2),SimpleBulkViscosity(0.1,15.0),HQdiffusion(0.2,1.5))
    dpt = pressure_derivative(1.0,Val(1),params.eos)
    @test Fluidum.viscosity(1.0,dpt,params.shear)!=0
    params=Fluidum.FluidProperties(eos,ZeroViscosity(),SimpleBulkViscosity(0.1,15.0),HQdiffusion(0.2,1.5))
    @test Fluidum.viscosity(1.0,dpt,params.shear)==0
    params=Fluidum.FluidProperties(eos,ZeroViscosity(),ZeroBulkViscosity(),HQdiffusion(0.2,1.5))
    @test Fluidum.bulk_viscosity(1.0,dpt,params.bulk)==0
    @test τ_bulk(1.0,dpt,params.bulk)==1
    params=Fluidum.FluidProperties(eos,ZeroViscosity(),ZeroBulkViscosity(),ZeroDiffusion())
    @test Fluidum.diffusion(1.0,dpt,params.diffusion)==0
end

#=
@testset "observables" begin 
    eos = Heavy_Quark()
    eos_HQ = HadronResonaceGas()
    params=Fluidum.FluidProperties(eos,QGPViscosity(0.2,0.2),SimpleBulkViscosity(0.1,15.0),HQdiffusion(0.2,1.5))
    fo = Fluidum.runFluidum_fo(Fluidum.Step_Intial_Condition(31.57,4),Fluidum.pQCD_Initial_Condition(1,70.,0.463),Fluidum.Trento_Intial_Condition(1),params,eos_HQ,0.4)
    obs = Fluidum.Observables(fo,0.200,params,0.156)
    Fluidum.save_observables(obs)
    @test isfile(Fluidum.get_filename(obs))
    rm(Fluidum.get_filename(obs))
    @test typeof(obs.yield_th)<:Float64
end

=#
#@testset "plots" begin
#    #Fluidum.plot_params(gui=true)
#    eos = Heavy_Quark()
#    obs=Fluidum.compute_observables(eos,1.5,Tfo=0.156,save = true)
#    println("ok...")
#    @test isfile(Fluidum.get_filename(obs))
#    println("ok...")
#    Fluidum.plot_spectra(obs,save=true)
#    println("ok...")
#    begin
#    result=Fluidum.runFluidum(eos,DsT=0,maxtime=5.)
#    Fluidum.plot_field(result,:temperature,tspan=(0.4,5.),save=true)
#    #Fluidum.plot_field(result,tspan=(0.4,5.),save=true)
#    println("ok...")
#    end
#    #rm(Fluidum.get_filename(obs)) 
#end

@testset "Cheb" begin
        function abs_matrix(diff, A)
           λ=eigen(A).values
           vec = eigen(A).vectors
        return vec*Diagonal(abs.(λ))*inv(vec)*diff
       end
       A = [1 1 ;1 1]
       A = A/max(eigen(A).values...)
       diff = [1.,1]
   
       B = copy(A)
       diff2 = copy(diff)
   
       Fluidum.cheb_flux!(diff,A,2)
     
       @test isapprox(diff,abs_matrix(diff2,B),rtol=0.2)
    
end

@testset "initialconditions" begin
    @test Fluidum.temperature(1,Fluidum.Trento_Intial_Condition(1))>0
    @test Fluidum.temperature(1,Fluidum.Trento_Intial_Condition(2))>Fluidum.temperature(1,Fluidum.Trento_Intial_Condition(1))
    @test Fluidum.dNdy(Fluidum.Step_Intial_Condition(31.57,4.0),Fluidum.pQCD_Initial_Condition(1,70.,0.463))==Fluidum.dNdy(Fluidum.Step_Intial_Condition(31.57,4.0),Fluidum.charm_pQCD())
    NQQ̄,err= quadgk(x->2*pi*x*0.4*thermodynamic(Fluidum.temperature(x,Fluidum.Trento_Intial_Condition(1)),Fluidum.fugacity(x,0.4,Fluidum.Step_Intial_Condition(31.57,4),Fluidum.pQCD_Initial_Condition(1,70.,0.463),Fluidum.Trento_Intial_Condition(1),Fluidum.HadronResonaceGas()),Fluidum.HadronResonaceGas()).pressure,0,30,rtol=0.00001)
    @test isapprox(NQQ̄,Fluidum.dNdy(Fluidum.Step_Intial_Condition(31.57,4.0),Fluidum.charm_pQCD())[1],rtol=0.2)
end

@testset "2d viscous " begin

fluidpropery=FluidProperties(FluiduMEoS(),SimpleShearViscosity(0.1,.1),ZeroBulkViscosity(),ZeroDiffusion())

# the convention here are T, ux, uy, piyy, pizz, pixy, piB this has to match with the matrix defined
twod_visc_hydro=Fields(
NDField((:ghost,:ghost),(:ghost,:ghost),:temperature),
NDField((:ghost,:ghost),(:ghost,:ghost),:ux),
NDField((:ghost,:ghost),(:ghost,:ghost),:uy),
NDField((:ghost,:ghost),(:ghost,:ghost),:piyy),
NDField((:ghost,:ghost),(:ghost,:ghost),:pizz),
NDField((:ghost,:ghost),(:ghost,:ghost),:pixy),
NDField((:ghost,:ghost),(:ghost,:ghost),:piB)
)

#we define a 2 cartesian grid form -25 to 25 50 point each dimension 
discretization=CartesianDiscretization(Fluidum.SymmetricInterval(50,25.),Fluidum.SymmetricInterval(50,25.))

# we prepare the field with the discretization
twod_visc_hydro_discrete=DiscreteFields(twod_visc_hydro,discretization,Float64)

 #we define some random intial condition 
function temperature(r)
       #0.4(1+0.3rand())/(exp(abs(r)-10)+1 )+0.01
       0.4/(exp(abs(r)-10)+1 )+0.01
end
#we set the array corresponding to the temperature 
phi=set_array((x,y)->temperature(hypot(x,y)),:temperature,twod_visc_hydro_discrete);


#this is the time span 
tspan=(0.4,20)


# this create the ODEProblem  and feed it diffenential equations 
res=oneshoot(twod_visc_hydro_discrete,Fluidum.matrix2d_visc!,fluidpropery,phi,tspan)

#plot the solution 

    result=res(20)
    @test all(isfinite.(result))
end 

#@testset "trento_initial_conditions" begin
#    IC=Fluidum.TrenTo_IC("Pb",0.0,1.5,0.5,0.5,1,0.0,6.4)
#    stats=Fluidum.TrenTo_IC_stats(10,1:50:100,20)
#    Fluidum.get_initial_profile(IC;statistics=stats)
#    @test isfile(Fluidum.assembleTrentoName(IC)[1])
#    @test isfile(Fluidum.assembleTrentoName(IC)[2])
#end

@testset "causality check" begin
    #define equation of state and transport coefficients
    ccbar = 30.
    #eos = Heavy_Quark(readresonancelist(), ccbar)
    eos = Heavy_Quark(readresonancelist(), ccbar)

    viscosity = QGPViscosity(0.1,0.2); #or, ZeroViscosity();
    bulk = SimpleBulkViscosity(0.083,15.0); #or, ZeroBulkViscosity();   

    diffusion_Ds = HQdiffusion(0.1, 1.5);
    params_Ds=Fluidum.FluidProperties(eos,viscosity,bulk,diffusion_Ds);
    #define a radial grid
    rmax = 20;
    gridpoints=500;
    grid_params = Fluidum.GridParameters(rmax, gridpoints); 

    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    #define the name and number of fields, with appropriate parity
    oned_visc_hydro = Fluidum.HQ_viscous_1d()
    disc_fields = DiscreteFields(oned_visc_hydro,disc,Float64) 
    #define a example function for the temperature
    function temperature(r)
           0.4*1/(exp(r/7)+1 )+0.0001
    end 

    #set the temperature field
    phi=set_array((x)->temperature(x),:temperature,disc_fields); 
    set_array!(phi,(x)->-3. *temperature(x),:α,disc_fields); 

    tspan = (0.4,10);
    field_results_Ds = Fluidum.oneshoot_debug(disc_fields, Fluidum.matrix1d_visc_HQ!, params_Ds, phi, tspan; reltol = 1e-8);
 
    iszero(field_results_Ds)

end

