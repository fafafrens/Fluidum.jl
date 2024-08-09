using Fluidum
using Test
using LinearAlgebra
using ForwardDiff
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

@testset "evolution" begin
    eos = Heavy_Quark()
    @test det(Fluidum.one_d_viscous_matrix([0.2,0.1,0,0,0,-0.1,0],2.,2.,0.5,0.5,0.5,0.5,0.5,0,0,0,1.,1.,0.1,0.1,0.1,1.,0)[1])!=0 
    @test typeof(Fluidum.runFluidum_fo(eos,Tfo=0.18,σ_temp=0.5).fo)<:FreezeOutResult
   # @test typeof(Fluidum.fields_evolution(eos,σ_temp=0.5,maxtime=2.0))<:Fields
end

@testset "observables" begin 
begin
    eos = Heavy_Quark()
    obs=Fluidum.compute_observables(eos,1.5,Tfo=0.18,σ_temp=0.5,save = true)
    @test isfile(Fluidum.get_filename(obs))
    @test typeof(obs.yield_th)<:Float64
    rm(Fluidum.get_filename(obs)) 
end
end

@testset "plots" begin
    begin
        Fluidum.plot_params(gui=true)
        eos = Heavy_Quark()
        obs=Fluidum.compute_observables(eos,1.5,Tfo=0.18,σ_temp=0.5,save = true)
        @test isfile(Fluidum.get_filename(obs))
        Fluidum.plot_spectra(obs,save=true)
    
    end
end

