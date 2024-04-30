using Fluidum
using Test

@testset "Fluidum.jl" begin
    # Write your tests here.
    @test IdealQCD()==IdealQCD{Float64}(3,2)
    @test isfinite(thermodynamic(1,IdealQCD()))
    @test isfinite(thermodynamic(0,IdealQCD()))
    @test isfinite(thermodynamic(10,IdealQCD(1,0)))

end
