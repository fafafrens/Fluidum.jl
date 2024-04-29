using Fluidum
using Test

@testset "Fluidum.jl" begin
    # Write your tests here.
    @test Fluidum.greet_your_package_name() == "Hello Fluidum!"
    @test Fluidum.greet_your_package_name() != "Hello world!"
end
