# The two-way medium + charm system (src/Matrix/1d_two_way_cilindrical.jl, 2d_two_way.jl).
using Fluidum, Test, StaticArrays, LinearAlgebra, OrdinaryDiffEq

# The closures are allocation-free on Julia >= 1.12. Older compilers allocate a few bytes per call;
# there the checks are recorded as skipped.
tw_noalloc(bytes) = VERSION >= v"1.12" ? (@test bytes == 0) : (@test_skip bytes == 0)

# P = a T^4 as a generic EOS (not a ConformalEOS), so the drivers take the general-EOS matrices
struct TWConformalAsGeneric <: Fluidum.EquationOfState
    a::Float64
end
Fluidum.thermodynamic(T, x::TWConformalAsGeneric) =
    Fluidum.Thermodynamic(x.a * T^4, (4x.a * T^3,), (12x.a * T^2,))

"Run a closure at (tau, X); return (A_t^-1 A_i..., A_t^-1 src)."
function tw_eval(m, ϕ, X, params; tau = 1.3)
    N = length(ϕ)
    A = ntuple(_ -> zeros(MMatrix{N,N,Float64}), length(X)); S = zeros(MVector{N,Float64})
    m(A, S, ϕ, tau, X, params)
    return map(SMatrix, A), SVector(S)
end

@testset "two-way system" begin
    F = Fluidum
    CH = charm_quark_list()
    CEOS = ConformalEOS(47.5)
    PARC = F.FluidProperties(CEOS, QGPViscosity(0.16, 0.2), ZeroBulkViscosity(), ConstDiffusion(0.116, 1.5))
    PARG = F.FluidProperties(TWConformalAsGeneric(F.conformal_aSB(CEOS)), QGPViscosity(0.16, 0.2), ZeroBulkViscosity(),
                             ConstDiffusion(0.116, 1.5))

    @testset "ConformalEOS" begin
        th = thermodynamic(0.3, CEOS); a = F.conformal_aSB(CEOS)
        @test th.pressure ≈ a * 0.3^4 && th.pressure_derivative[1] ≈ 4a * 0.3^3 && th.pressure_hessian[1] ≈ 12a * 0.3^2
        @test thermodynamic(1, ConformalEOS(40)).pressure isa Float64
    end

    @testset "conformal matrices ≡ general-EOS matrices on P = a T^4" begin
        # the two exports of the same derivation; the bulk field of the general layout is held at 0
        k1 = [1, 2, 3, 4, 6, 7, 8, 9, 10]
        ϕ = SVector(0.3, 0.2, 0.001, -0.002, -2.0, 1e-4, 1e-4, -5e-5, 2e-5)
        ϕg = SVector(ϕ[1:4]..., 0.0, ϕ[5:9]...)
        for edge in (nothing, ColdEdge())
            (Ac,), Sc = tw_eval(matrix1d_two_way!(; species = CH, edge), ϕ, (2.7,), PARC)
            (Ag,), Sg = tw_eval(matrix1d_two_way!(; species = CH, edge), ϕg, (2.7,), PARG)
            @test maximum(abs.(Ag[k1, k1] .- Ac)) <= 1e-12 * maximum(abs, Ac)
            @test maximum(abs.(Sg[k1] .- Sc)) <= 1e-12 * maximum(abs, Sc)
        end
        k2 = [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14]
        ϕ2 = SVector(0.3, 0.2, -0.1, 0.001, -0.002, 0.0003, -2.0, 0.01, -0.005, 0.001, -0.0005, 0.0002, 0.0001)
        ϕ2g = SVector(ϕ2[1:6]..., 0.0, ϕ2[7:13]...)
        (Acx, Acy), Sc2 = tw_eval(matrix2d_two_way!(; species = CH), ϕ2, (3.0, 1.0), PARC)
        (Agx, Agy), Sg2 = tw_eval(matrix2d_two_way!(; species = CH), ϕ2g, (3.0, 1.0), PARG)
        @test maximum(abs.(Agx[k2, k2] .- Acx)) <= 1e-12 * maximum(abs, Acx)
        @test maximum(abs.(Agy[k2, k2] .- Acy)) <= 1e-12 * maximum(abs, Acy)
        @test maximum(abs.(Sg2[k2] .- Sc2)) <= 1e-12 * maximum(abs, Sc2)
    end

    @testset "the second moment does not enter the medium and first-moment rows" begin
        a = F.conformal_aSB(CEOS)
        ϕ = SVector(0.3, 0.2, 0.001, -0.002, -2.0, 1e-4, 1e-4, -5e-5, 2e-5)
        c = F.charm_eos(0.3, -2.0, CH)
        At, Ar, src = F.one_d_two_way_matrix_conformal(ϕ, 1.3, 2.7, a, 0.1, 0.5, 0.6, values(c)..., 0.4)
        @test all(iszero, At[1:6, 7:9]) && all(iszero, Ar[1:6, 7:9])
        ϕb = SVector(ϕ[1:6]..., 3e-4, 2e-4, -1e-4)
        @test F.one_d_two_way_matrix_conformal(ϕb, 1.3, 2.7, a, 0.1, 0.5, 0.6, values(c)..., 0.4)[3][1:6] == src[1:6]
    end

    @testset "Bjorken: the charm number is conserved" begin
        # uniform state, u^r = 0: dφ/dτ = -A_t^-1 src, and d(τ n)/dτ = 0
        for (T, α) in ((0.3, -2.0), (0.2, 1.0))
            ϕ = SVector(T, 0.0, 0.002, -0.004, α, 0.0, 0.0, 0.0, 0.0)
            _, S = tw_eval(matrix1d_two_way!(; species = CH, dmn_eps = 0.0), ϕ, (2.7,), PARC; tau = 1.3)
            c = F.charm_eos(T, α, CH)
            @test 1.3 * (c.dndT * (-S[1]) + c.dnda * (-S[5])) + c.n ≈ 0 atol = 1e-12 * c.n
        end
    end

    @testset "cold edge" begin
        e = ColdEdge()
        @test F.edge_factor(0.2, e) == 1.0 && F.edge_factor(0.156, e) == 1.0
        @test F.edge_factor(0.12, e) == e.floor && 0 < F.edge_factor(0.146, e) < 1
        @test F.edge_factor(0.12, nothing) == 1.0
        @test_throws ErrorException ColdEdge(; floor = 0.0)
        ϕ = SVector(0.3, 0.2, 0.001, -0.002, -2.0, 1e-4, 1e-4, -5e-5, 2e-5)        # above T_edge: untouched
        @test tw_eval(matrix1d_two_way!(; species = CH, edge = nothing), ϕ, (2.7,), PARC) ==
              tw_eval(matrix1d_two_way!(; species = CH), ϕ, (2.7,), PARC)
        ϕc = SVector(0.13, 0.2, 0.001, -0.002, -2.0, 1e-4, 1e-4, -5e-5, 2e-5)       # below: D_s is reduced
        @test tw_eval(matrix1d_two_way!(; species = CH, edge = nothing), ϕc, (2.7,), PARC) !=
              tw_eval(matrix1d_two_way!(; species = CH), ϕc, (2.7,), PARC)
    end

    @testset "closures: zero allocations" begin
        allocs(m, ϕ, X, P) = (tw_eval(m, ϕ, X, P); N = length(ϕ);
            A = ntuple(_ -> zeros(MMatrix{N,N,Float64}), length(X)); S = zeros(MVector{N,Float64});
            m(A, S, ϕ, 1.0, X, P); @allocated m(A, S, ϕ, 1.0, X, P))
        ϕ = SVector(0.3, 0.2, 0.001, -0.002, -2.0, 0.01, 0.001, -0.0005, 0.0002)
        tw_noalloc(allocs(matrix1d_two_way!(; species = CH), ϕ, (3.0,), PARC))
        ϕ2 = SVector(0.3, 0.2, -0.1, 0.001, -0.002, 0.0003, -2.0, 0.01, -0.005, 0.001, -0.0005, 0.0002, 0.0001)
        tw_noalloc(allocs(matrix2d_two_way!(; species = CH), ϕ2, (3.0, 1.0), PARC))
    end

    @testset "two_way_problem: the corona update, and a short 1+1D run" begin
        disc = CartesianDiscretization(F.OriginInterval(60, 12.0))
        d = DiscreteFields(two_way_conformal_1d(), disc, Float64)
        Tr(r) = 0.12 + 0.25 / (1 + exp((r - 5) / 0.6))
        phi = set_array(Tr, :temperature, d)
        # one charm quark: n ∝ e^α
        set_array!(phi, r -> log((0.2 * exp(-(r / 3)^2) + 1e-12) / F.charm_eos(Tr(r), 0.0, CH).n), :α, d)
        p1 = two_way_problem(d, PARC, phi, (0.4, 3.0); species = CH)
        p0 = two_way_problem(d, PARC, phi, (0.4, 3.0); species = CH, edge = ColdEdge(; corona = false))
        du1 = similar(phi); du0 = similar(phi); p1.f(du1, phi, nothing, 0.4); p0.f(du0, phi, nothing, 0.4)
        hot = [phi[1, i] >= ColdEdge().T_edge for i in axes(phi, 2)]
        # above T_edge nothing changes; below it only the α rate (row 5) does
        @test du1[:, hot] == du0[:, hot] && du1[[1, 2, 3, 4, 6, 7, 8, 9], :] == du0[[1, 2, 3, 4, 6, 7, 8, 9], :]
        @test du1[5, .!hot] != du0[5, .!hot]
        rr = [disc.grid[i][1] for i in eachindex(disc.grid)]
        allocs(f, a...) = (f(a...); @allocated f(a...))        # in a function: no global-scope boxing
        tw_noalloc(allocs(F.corona_update_1d!, du1, phi, 0.4, rr, 0.2, 5, CH, ColdEdge(), 1e-8, 200.0))
        sol = solve(p1, Tsit5(); reltol = 1e-8, abstol = 1e-10, dtmax = 0.01)
        @test string(sol.retcode) == "Success" && all(isfinite, sol.u[end])
    end

    @testset "field layouts" begin
        @test length(F.names(two_way_1d())) == 10 && length(F.names(two_way_conformal_1d())) == 9
        @test length(F.names(two_way_2d())) == 14 && length(F.names(two_way_conformal_2d())) == 13
    end
end
