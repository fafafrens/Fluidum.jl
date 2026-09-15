using Test
using StaticArrays

include("thermodynamics_v2.jl")
using .FluidumThermodynamicsPrototype

@testset "Thermodynamics v2 prototype" begin
    χ0 = @SMatrix [2.0 0.3; 0.3 1.5]
    eos = PolynomialMultiChargeEOS(3.0, χ0)

    x = GrandCanonicalPoint(0.2, 0.03, -0.01)
    th = thermodynamic(eos,x)

    @test length(th.gradient) == 3
    @test size(th.hessian) == (3,3)
    @test charge_densities(th) ≈ SVector(th.gradient[2], th.gradient[3])
    @test susceptibilities(th) ≈ th.hessian[2:3,2:3]

    ε = energy_density(x,th)
    @test isfinite(ε)

    can = CanonicalPoint(x.T, charge_densities(th))
    xcan = grandcanonical(eos, can, @SVector [0.0, 0.0])
    @test xcan.T ≈ x.T
    @test xcan.μ ≈ x.μ atol=1e-10 rtol=1e-10

    con = ConservedPoint(entropy_density(th), charge_densities(th))
    xfull = invert(eos, con, GrandCanonicalPoint(0.18, 0.0, 0.0))
    @test xfull.T ≈ x.T atol=1e-10 rtol=1e-10
    @test xfull.μ ≈ x.μ atol=1e-10 rtol=1e-10

    κ = @SMatrix [0.4 0.05; 0.05 0.3]
    τ = @SMatrix [1.0 0.0; 0.0 1.2]
    tm = TransportModel(
        ConstantEtaOverS(0.16, 5.0),
        ConstantZetaOverS(0.03, 10.0),
        ConstantChargeTransport(κ,τ),
    )
    fluid = FluidModel(eos,tm)

    th2,tr = transport(fluid,x)
    @test th2.pressure ≈ th.pressure
    @test tr.shear.eta > 0
    @test tr.bulk.zeta > 0
    @test tr.charge.kappa == κ
    @test tr.charge.tau == τ

    # EOS algebra lives at EOS level, not ThermodynamicState level.
    eos2 = 2.0*eos + (-eos)
    thalg = thermodynamic(eos2,x)
    @test thalg.pressure ≈ th.pressure
    @test thalg.gradient ≈ th.gradient
    @test thalg.hessian ≈ th.hessian
end
