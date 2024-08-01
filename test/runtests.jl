using Fluidum
using Test
using LinearAlgebra

@testset "Fluidum.jl" begin
    # Write your tests here.
    @test IdealQCD()==IdealQCD(3,2)
    @test isfinite(thermodynamic(1,IdealQCD()).pressure[1])
    @test isfinite(thermodynamic(0,IdealQCD()).pressure_derivative[1])
  
    @test length(thermodynamic(10,IdealQCD(1,0)).pressure)==1
    @test length(thermodynamic(10,IdealQCD(1,0)).pressure_derivative)==1
    @test length(thermodynamic(10,IdealQCD(1,0)).pressure_hessian)==1

    @test isfinite(free_charm(1,0,Heavy_Quark())[1])
    @test pressure(1,Heavy_Quark())==pressure(1,FluiduMEoS())
    @test pressure(1,Heavy_Quark())==thermodynamic(1,FluiduMEoS()).pressure[1]

    @test det(Fluidum.one_d_viscous_matrix([0.2,0.1,0,0,0,-0.1,0],2.,2.,0.5,0.5,0.5,0.5,0.5,0,0,0,1.,1.,0.1,0.1,0.1,1.,0)[1])!=0
    
    @test Fluidum.fluidproperties(FluiduMEoS()).eos==FluiduMEoS()
end


