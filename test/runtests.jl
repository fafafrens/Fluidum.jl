using Fluidum
using Test
using LinearAlgebra
using ForwardDiff
using QuadGK
#=
@testset "equation of state" begin
    @test IdealQCD()==IdealQCD(3,2)
    @test isfinite(thermodynamic(1,IdealQCD()).pressure[1])
    @test isfinite(thermodynamic(0,IdealQCD()).pressure_derivative[1])
  
    @test length(thermodynamic(10,IdealQCD(1,0)).pressure)==1
    @test length(thermodynamic(10,IdealQCD(1,0)).pressure_derivative)==1
    @test length(thermodynamic(10,IdealQCD(1,0)).pressure_hessian)==1

    @test isfinite(free_charm(1,0,Heavy_Quark())[1])
    @test round(pressure(1,Heavy_Quark()), sigdigits=10)==round(pressure(1,FluiduMEoS()), sigdigits=10)
    @test round(pressure(1,Heavy_Quark()), sigdigits=10)==round(thermodynamic(1,FluiduMEoS()).pressure[1], sigdigits=10)
    @test pressure(1,Heavy_Quark()) ≈ pressure(1,FluiduMEoS()) atol=0.01   
end
=#
@testset "fluid_properties" begin
    eos = FluiduMEoS()
    @test det(Fluidum.one_d_viscous_matrix([0.2,0.1,0,0,0,-0.1,0],2.,2.,0.5,0.5,0,0,0,1.,1.,0.1,0.1,0.1,1.,0)[1])!=0 
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
#=
@testset "plots" begin
    #Fluidum.plot_params(gui=true)
    eos = Heavy_Quark()
    obs=Fluidum.compute_observables(eos,1.5,Tfo=0.156,save = true)
    println("ok...")
    @test isfile(Fluidum.get_filename(obs))
    println("ok...")
    Fluidum.plot_spectra(obs,save=true)
    println("ok...")
    begin
    result=Fluidum.runFluidum(eos,DsT=0,maxtime=5.)
    Fluidum.plot_field(result,:temperature,tspan=(0.4,5.),save=true)
    #Fluidum.plot_field(result,tspan=(0.4,5.),save=true)
    println("ok...")
    end
    #rm(Fluidum.get_filename(obs)) 
end
=#
#=@testset "Cheb" begin
    function abs_matrix(diff, A)
        λ=eigen(A).values
        vec = eigen(A).vectors
        return vec*Diagonal(abs.(λ))*inv(vec)*diff
    end
    A = rand(2,2)
    A = A/max(eigen(A).values...)
    diff = rand(2)

    B = copy(A)
    diff2 = copy(diff)

    Fluidum.cheb_flux!(diff,A,15)
    @test isapprox(diff,abs_matrix(diff2,B),rtol=0.2)
    
end=# #FC 21.10.24: cheb test not always passing

@testset "initialconditions" begin
    @test Fluidum.temperature(1,Fluidum.Trento_Intial_Condition(1))>0
    @test Fluidum.temperature(1,Fluidum.Trento_Intial_Condition(2))>Fluidum.temperature(1,Fluidum.Trento_Intial_Condition(1))
    @test Fluidum.dNdy(Fluidum.Step_Intial_Condition(31.57,4.0),Fluidum.pQCD_Initial_Condition(1,70.,0.463))==Fluidum.dNdy(Fluidum.Step_Intial_Condition(31.57,4.0),Fluidum.charm_pQCD())
    NQQ̄,err= quadgk(x->2*pi*x*0.4*thermodynamic(Fluidum.temperature(x,Fluidum.Trento_Intial_Condition(1)),Fluidum.fugacity(x,0.4,Fluidum.Step_Intial_Condition(31.57,4),Fluidum.pQCD_Initial_Condition(1,70.,0.463),Fluidum.Trento_Intial_Condition(1),Fluidum.HadronResonaceGas()),Fluidum.HadronResonaceGas()).pressure,0,30,rtol=0.00001)
    @test isapprox(NQQ̄,Fluidum.dNdy(Fluidum.Step_Intial_Condition(31.57,4.0),Fluidum.charm_pQCD())[1],rtol=0.2)
end
#end


