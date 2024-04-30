using Fluidum
using Test

@testset "Fluidum.jl" begin
    # Write your tests here.
    @test IdealQCD()==IdealQCD(3,2)
    @test isfinite(first(thermodynamic(1,IdealQCD())))
    @test isfinite(thermodynamic(0,IdealQCD())[2])
    @test length(thermodynamic(10,IdealQCD(1,0)))==3
    @test length(first(thermodynamic(10,IdealQCD(1,0))))==1
    @test length(last(thermodynamic(10,IdealQCD(1,0))))==1

end
