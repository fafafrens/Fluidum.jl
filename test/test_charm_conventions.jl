# Conventions of the charm list, the linear D_s T law and the dict initializer

const besseli = Fluidum.besseli

@testset "charm conventions" begin

    @testset "canonical factor by list type" begin
        q = charm_quark_list().particle_list                     # one open-charm state, q = 1
        can = Fluidum.HadronResonaceGas_ccbar(q, 30.0)
        gc  = HadronResonaceGas_ccbar_grandcanonical(q, 30.0)
        f = besseli(1, 15.0) / besseli(0, 15.0)
        @test Fluidum.canonical_factor(can) ≈ f
        @test Fluidum.canonical_factor(gc) == 1
        # the density and its derivatives carry the factor for the canonical list only
        for (T, α) in ((0.16, 0.0), (0.25, 1.3))
            tc, tg = Fluidum.thermodynamic(T, α, can), Fluidum.thermodynamic(T, α, gc)
            @test tc.pressure ≈ f * tg.pressure
            @test all(tc.pressure_derivative .≈ f .* tg.pressure_derivative)
        end
    end

    @testset "linear D_s T law frozen below Tc" begin
        a, b = 1.765, -0.159
        d = LinearDiffusion(a, b, 1.5)
        @test d.Tc == 0.180
        @test Fluidum.DsT(d, 0.10) ≈ a * 0.180 + b
        @test Fluidum.DsT(d, 0.30) ≈ a * 0.30 + b
        d156 = LinearDiffusion(a, b, 1.5; Tc = 0.156)
        @test Fluidum.DsT(d156, 0.10) ≈ a * 0.156 + b
        @test Fluidum.DsT(d156, 0.17) ≈ a * 0.17 + b
        @test Fluidum.DsT(d156, 0.10; Tc = 0.2) ≈ a * 0.2 + b       # the keyword still wins
    end

    @testset "dict initializer on a prescribed background" begin
        grid = (rmax = 10.0, gridpoints = 40)
        f = Dict{Symbol,Function}(:α => r -> 0.3 - 0.01r, :nur => r -> 1e-3r,
                                  :piphiphi_Q => r -> 2e-4, :pietaeta_Q => r -> -4e-4, :piB_Q => r -> -1e-4)
        init = Fluidum.initialize_fields(f, Fluidum.HQ_const_BG_second_moment_init, grid)
        # the same array built by hand, field by field
        disc = Fluidum.CartesianDiscretization(OriginInterval(grid.gridpoints, grid.rmax))
        df = Fluidum.DiscreteFields(Fluidum.HQ_const_BG_second_moment_init(), disc, Float64)
        ref = Fluidum.set_array(f[:α], :α, df)
        Fluidum.set_array!(ref, f[:nur], :nur, df)
        for k in (:piphiphi_Q, :pietaeta_Q, :piB_Q)
            Fluidum.set_array!(ref, f[k], k, df)
        end
        @test init.initial_field == ref
    end
end
