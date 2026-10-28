# Charm on a prescribed background: the charm EOS inputs, the first and second moment, the 2-D -> 1-D
# reduction, zero allocations.
using Fluidum, Test, StaticArrays, LinearAlgebra, QuadGK

const HQM_CHARM_FILE = let p = joinpath(mktempdir(), "charm.data")
    write(p, "# Name Mass Gamma deg Spin Isospin I3 Nq Ns Naq Nas Nc Nac MC\n" *
             "Charm 1.5 0.0 6 0.5 0. 0. 0. 0. 0. 0. 1. 0. 4\n")
    p
end
# one charm quark, grand canonical (no canonical factor), as the charm matrices assume
const HQM_EOS = Fluidum.Heavy_Quark(hadron_list = HadronResonaceGas_ccbar_grandcanonical(
                    readresonancelist(name_file = HQM_CHARM_FILE), 1.0))
const HQM_PARAMS = Fluidum.FluidProperties(HQM_EOS, QGPViscosity(0.1, 0.2),
                                           SimpleBulkViscosity(0.083, 15.0), ConstDiffusion(0.11634, 1.5))
const HQM_BG1 = (tau, r) -> (0.31, 0.42, -0.055, -0.021, 0.13, 0.07)
const HQM_BG2 = (tau, x, y) -> (0.3, 0.4, -0.2, -0.05, -0.02, 0.01, 0.05, 0.1, 0.02, -0.03, 0.04, 0.08)
const HQM_PHI5 = SVector(-2.0, 0.031, 1.0e-3, -4.0e-4, 2.0e-4)
const HQM_PHI7 = SVector(-2.0, 0.01, -0.02, 1e-3, 2e-4, -5e-4, 3e-4)

"Run a closure at (tau, X); return (A_t^-1 A_i..., A_t^-1 src)."
function hqm_eval(m, ϕ, X; tau = 2.3, params = HQM_PARAMS)
    N = length(ϕ)
    A = ntuple(_ -> zeros(MMatrix{N,N,Float64}), length(X)); S = zeros(MVector{N,Float64})
    m(A, S, ϕ, tau, X, params)
    return map(SMatrix, A), SVector(S)
end
# The drivers are allocation-free on Julia >= 1.12. Older compilers allocate a few bytes per call in
# these closures; there the checks are recorded as skipped.
hqm_noalloc(bytes) = VERSION >= v"1.12" ? (@test bytes == 0) : (@test_skip bytes == 0)
function hqm_alloc(m, ϕ, X)
    N = length(ϕ)
    A = ntuple(_ -> zeros(MMatrix{N,N,Float64}), length(X)); S = zeros(MVector{N,Float64})
    m(A, S, ϕ, 1.2, X, HQM_PARAMS); m(A, S, ϕ, 1.2, X, HQM_PARAMS)
    return @allocated m(A, S, ϕ, 1.2, X, HQM_PARAMS)
end

# lifts a 1-D background bg1d(tau, r) -> (T, ur, dtT, drT, drur, dtur) to the 2-D one
function hqm_axisymmetric(bg1d)
    return (tau, x, y) -> begin
        r = sqrt(x^2 + y^2); cx, cy = x / r, y / r
        T, ur, dtT, drT, drur, dtur = bg1d(tau, r)
        s = drur - ur / r      # d_j u^i = (u^r' - u^r/r) x^i x_j / r² + (u^r/r) δ^i_j
        (T, ur * cx, ur * cy, dtT, drT * cx, drT * cy, dtur * cx, s * cx * cx + ur / r, s * cx * cy,
         dtur * cy, s * cy * cx, s * cy * cy + ur / r)
    end
end

@testset "HQ on a prescribed background" begin

@testset "charm_eos" begin
    l = charm_quark_list()
    # the canonical factor: in HadronResonaceGas_ccbar, not in the grand-canonical charm_quark_list
    lc = HadronResonaceGas_ccbar(l.particle_list, 30.0)
    @test thermodynamic(0.3, -2.0, lc).pressure ≈
          Fluidum.besseli(1, 15.0) / Fluidum.besseli(0, 15.0) * thermodynamic(0.3, -2.0, l).pressure rtol = 1e-14
    for (T, α) in ((0.3, -2.0), (0.16, 0.5), (0.5, -10.0))
        c = Fluidum.charm_eos(T, α, l)
        z = 1.5 / T
        @test c.h ≈ 1.5 * Fluidum.besselk(3, z) / Fluidum.besselk(2, z) rtol = 1e-13
        @test c.cQ ≈ 1.5 * Fluidum.besselk(4, z) / Fluidum.besselk(3, z) rtol = 1e-13
        # tau_n = D_s h/T, the relaxation time of the Fick closure, for one quark
        @test Fluidum.τ_diffusion_hadron(T, α, HQM_EOS, HQM_PARAMS.diffusion) ≈
              Fluidum.DsT(HQM_PARAMS.diffusion, T) / T / Fluidum.fmGeV * c.h / T rtol = 1e-12
        δ = 1e-5
        @test c.dhdT ≈ (Fluidum.charm_eos(T + δ, α, l).h - Fluidum.charm_eos(T - δ, α, l).h) / 2δ rtol = 1e-7
        @test c.dcQdT ≈ (Fluidum.charm_eos(T + δ, α, l).cQ - Fluidum.charm_eos(T - δ, α, l).cQ) / 2δ rtol = 1e-7
    end
    # three unit-charge species: the species sums against a momentum quadrature
    L3 = HadronResonaceGas_ccbar_grandcanonical(Fluidum.StructArray([
             Fluidum.Particle("D", 1.8646, 0.0, 4, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 421),
             Fluidum.Particle("Dst", 2.007, 0.0, 12, 1.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 423),
             Fluidum.Particle("Lc", 2.2865, 0.0, 4, 0.5, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 4122)]), 1.0)
    # I_nq by quadrature; ∂_T I_nq|_α = I_{n+1,q}/T² exactly
    function Inq(n, q, T, α, m, g)
        z = m / T; tm = acosh(1 + 900 / z)
        f(t) = m * cosh(t) * (m * sinh(t))^(2 + 2q) * (m * cosh(t))^(n - 2q - 1) * exp(-z * (cosh(t) - 1))
        return g * exp(α - z) / (2π^2 * prod(1:2:(2q+1))) * quadgk(f, 0.0, tm / 4, tm / 2, tm; rtol = 1e-14)[1]
    end
    worst = 0.0
    for T in (0.16, 0.25, 0.4), α in (-3.0, 1.0)
        S(n, q) = sum(Inq(n, q, T, α, p.Mass, p.Degeneracy) for p in L3.particle_list)
        h = S(3, 1) / S(2, 1); dh = (S(4, 1) - h * S(3, 1)) / T^2 / S(2, 1)
        C = S(5, 2) / S(4, 2); dC = (S(6, 2) - C * S(5, 2)) / T^2 / S(4, 2)
        A = S(4, 1) / S(2, 1)
        c = Fluidum.charm_eos(T, α, L3)
        worst = max(worst, abs(c.h / h - 1), abs(c.dhdT / dh - 1), abs(c.cQ / C - 1), abs(c.dcQdT / dC - 1),
                    abs(c.aQ / A - 1), abs(T * c.dI31dT / c.I31 / (T * S(4, 1) / (T^2 * S(3, 1))) - 1))
    end
    @test worst < 1e-8
    hqm_noalloc((Fluidum.charm_eos(0.3, -1.0, L3); @allocated Fluidum.charm_eos(0.3, -1.0, L3)))
end

@testset "second moment: rows 1-2 are the consistent first moment, and m2 is passive" begin
    (A2,), S2 = hqm_eval(matrix1d_visc_HQ_BG_consistent!(; background_fields = HQM_BG1), HQM_PHI5[1:2], (3.1,))
    (A5,), S5 = hqm_eval(matrix1d_visc_HQ_BG_consistent_m2!(; background_fields = HQM_BG1), HQM_PHI5, (3.1,))
    @test isapprox(A5[1:2, 1:2], A2; rtol = 1e-12) && isapprox(S5[1:2], S2; rtol = 1e-12)
    # same (α, ν), other moments: rows 1-2 do not change (to roundoff, after a dense A_t^-1)
    ϕb = SVector(HQM_PHI5[1], HQM_PHI5[2], -3.0e-3, 7.0e-4, 1.1e-3)
    (Ab,), Sb = hqm_eval(matrix1d_visc_HQ_BG_consistent_m2!(; background_fields = HQM_BG1), ϕb, (3.1,))
    @test maximum(abs, A5[1:2, 3:5]) <= 1e-14 * maximum(abs, A5)
    @test isapprox(S5[1:2], Sb[1:2]; rtol = 1e-13) && isapprox(A5[1:2, 1:2], Ab[1:2, 1:2]; rtol = 1e-13)
    # 2-D
    (Bx, By), T3 = hqm_eval(matrix2d_HQ_BG_consistent!(; background_fields = HQM_BG2), HQM_PHI7[1:3], (3.0, -1.0))
    (Cx, Cy), T7 = hqm_eval(matrix2d_HQ_BG_m2!(; background_fields = HQM_BG2), HQM_PHI7, (3.0, -1.0))
    @test isapprox(Cx[1:3, 1:3], Bx; rtol = 1e-12) && isapprox(Cy[1:3, 1:3], By; rtol = 1e-12)
    @test isapprox(T7[1:3], T3; rtol = 1e-12)
    @test maximum(abs, Cx[1:3, 4:7]) <= 1e-14 * maximum(abs, Cx)
end

@testset "2-D first moment reduces to the 1-D one (axisymmetric, on the x axis)" begin
    bg2 = hqm_axisymmetric(HQM_BG1)
    for (m1, m2) in ((matrix1d_visc_HQ_BG!(; background_fields = HQM_BG1), matrix2d_HQ_BG!(; background_fields = bg2)),
                     (matrix1d_visc_HQ_BG_consistent!(; background_fields = HQM_BG1),
                      matrix2d_HQ_BG_consistent!(; background_fields = bg2)))
        (Ar,), Sr = hqm_eval(m1, SVector(-2.0, 0.031), (3.1,))
        (Ax, Ay), Sx = hqm_eval(m2, SVector(-2.0, 0.031, 0.0), (3.1, 0.0))
        @test isapprox(Ax[1:2, 1:2], Ar; rtol = 1e-12)
        # the cylindrical nu^r/r is d_y nu^y = nu/r on the axis, so it enters through A_y
        S2 = Sx + Ay * SVector(0.0, 0.0, 0.031 / 3.1)
        @test isapprox(S2[1:2], Sr; rtol = 1e-12) && abs(S2[3]) <= 1e-14 * maximum(abs, Sr)
    end
end

@testset "zero allocations per cell" begin
    hqm_noalloc(hqm_alloc(matrix1d_visc_HQ_BG!(; background_fields = HQM_BG1), SVector(-2.0, 0.01), (3.0,)))
    hqm_noalloc(hqm_alloc(matrix1d_visc_HQ_BG_consistent!(; background_fields = HQM_BG1), SVector(-2.0, 0.01), (3.0,)))
    hqm_noalloc(hqm_alloc(matrix1d_visc_HQ_BG_consistent_m2!(; background_fields = HQM_BG1), HQM_PHI5, (3.0,)))
    hqm_noalloc(hqm_alloc(matrix2d_HQ_BG!(; background_fields = HQM_BG2), HQM_PHI7[1:3], (3.0, -1.0)))
    hqm_noalloc(hqm_alloc(matrix2d_HQ_BG_consistent!(; background_fields = HQM_BG2), HQM_PHI7[1:3], (3.0, -1.0)))
    hqm_noalloc(hqm_alloc(matrix2d_HQ_BG_m2!(; background_fields = HQM_BG2), HQM_PHI7, (3.0, -1.0)))
end

@testset "field layouts" begin
    @test length(Fluidum.names(HQ_const_BG_init())) == 2 && length(Fluidum.names(HQ_const_BG_second_moment_init())) == 5
    @test length(Fluidum.names(HQ_2p1d_BG_init())) == 3 && length(Fluidum.names(HQ_2p1d_BG_m2_init())) == 7
    @test length(Fluidum.names(HQ_viscous_2d())) == 10
end

end
