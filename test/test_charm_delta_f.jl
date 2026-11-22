# The off-equilibrium correction of the charm Cooper–Frye spectrum: CharmDeltaF, spectra_internal(...; delta_f)
# (src/fluidevo/spectra_fastreso_dict_HQ.jl), and the δf of the integrated spectrum (spectra_fastreso2_dict_HQ.jl).
using Fluidum, Test, StaticArrays, QuadGK, Random

const DF_CH = charm_quark_list()           # m = 1.5 GeV, g = 6
const DF_M = 1.5
const DF_DEG = 6
const DF_FMGEV = 5.068                     # the constant internal_thermal_spectra uses

# f₀ δ integrated over η and φ numerically, δ(k_l, k_φ, k_η) given; normalisation of internal_thermal_spectra
function df_brute(δ, pt, r, t, dra, dta, ur, T, α; nη = 4000, nφ = 512)
    mt = sqrt(DF_M^2 + pt^2); ut = sqrt(1 + ur^2)
    ηmax = acosh(1 + 60T / (ut * mt)); dη = 2ηmax / nη; dφ = 2π / nφ
    s = 0.0
    for i in 0:nη, j in 0:nφ-1
        η = -ηmax + i * dη; w = (i == 0 || i == nη) ? 0.5 : 1.0
        φ = j * dφ; E = mt * cosh(η); cφ = cos(φ)
        kl = ut * pt * cφ - ur * E; kφ = pt * sin(φ); kη = mt * sinh(η)
        s += w * (dta * pt * cφ - dra * E) * exp(α - (ut * E - ur * pt * cφ) / T) * δ(kl, kφ, kη)
    end
    return s * dη * dφ * t * r / (2π)^3 * DF_FMGEV^3 * DF_DEG
end

# a random node: (T, ur, α), nur of order 0.3 n, and piphiphi_Q, pietaeta_Q, piB_Q of order 0.1 n T
function df_node(rng)
    T = 0.15 + 0.15rand(rng); ur = 1.5rand(rng); α = -3 + 3rand(rng)
    n = Fluidum.charm_eos(T, α, DF_CH).n
    piphiphi_Q, pietaeta_Q, piB_Q = (0.2n * T) .* (2 .* rand(rng, 3) .- 1)
    return (; T, ur, α, n, νr = 0.3n * (2rand(rng) - 1), piphiphi_Q, pietaeta_Q, piB_Q,
            r = 1 + 6rand(rng), t = 1 + 9rand(rng), dra = 2rand(rng) - 1, dta = rand(rng))
end
df_f0(s, pt; ν = 0.0) = Fluidum.internal_thermal_spectra(pt, DF_M, s.r, s.t, s.dra, s.dta, s.ur, s.T, s.α, ν; deg = DF_DEG)
df_scale(s, pt) = Fluidum.internal_thermal_spectra(pt, DF_M, s.r, s.t, -abs(s.dra), abs(s.dta), s.ur, s.T, s.α, 0.0;
                                                   deg = DF_DEG)
df_d2(s, pt) = Fluidum.internal_df2_spectra(pt, DF_M, s.r, s.t, s.dra, s.dta, s.ur, s.T, s.α, s.piphiphi_Q, s.pietaeta_Q, s.piB_Q, DF_CH;
                                            deg = DF_DEG)

@testset "δf in the charm spectrum" begin
    F = Fluidum

    @testset "closed forms ≡ η×φ quadrature" begin
        rng = MersenneTwister(20261027); w1 = w2 = 0.0
        for _ in 1:12
            s = df_node(rng); pt = 0.1 + 2.9rand(rng); ut = sqrt(1 + s.ur^2)
            # δf₁: ν_l k_l/(nT) with ν_l = ν^r/u^τ
            d1 = df_f0(s, pt; ν = s.νr / s.n) - df_f0(s, pt)
            b1 = df_brute((kl, kφ, kη) -> s.νr / ut * kl / (s.n * s.T), pt, s.r, s.t, s.dra, s.dta, s.ur, s.T, s.α)
            w1 = max(w1, abs(d1 - b1) / df_scale(s, pt))
            c = F.df2_coefficients(s.T, s.α, s.piphiphi_Q, s.pietaeta_Q, s.piB_Q, DF_CH)
            b2 = df_brute((kl, kφ, kη) -> c.Ar * kl^2 + c.Aφ * kφ^2 + c.Aη * kη^2 + c.d,
                          pt, s.r, s.t, s.dra, s.dta, s.ur, s.T, s.α)
            w2 = max(w2, abs(df_d2(s, pt) - b2) / df_scale(s, pt))
        end
        @test w1 < 1e-7
        @test w2 < 1e-7
    end

    @testset "δf₂ matching in the rest frame: the charm stress reproduced, n unchanged" begin
        for T in (0.156, 0.25, 0.4), α in (-2.0, 0.5)
            n = F.charm_eos(T, α, DF_CH).n
            piphiphi_Q, pietaeta_Q, piB_Q = -0.05n * T, -0.08n * T, 0.07n * T
            c = F.df2_coefficients(T, α, piphiphi_Q, pietaeta_Q, piB_Q, DF_CH)
            f0(k) = DF_DEG * exp(α - sqrt(DF_M^2 + k^2) / T) * F.fmGeV3 * k^2 / (2π^2)
            E(k) = sqrt(DF_M^2 + k^2)
            Σ = c.Ar + c.Aφ + c.Aη
            # angular averages: ⟨n_i²⟩ = 1/3, ⟨n_i⁴⟩ = 1/5, ⟨n_i² n_j²⟩ = 1/15
            q(g) = quadgk(g, 0, 40T + 10; rtol = 1e-13)[1]
            @test q(f0) ≈ n rtol = 1e-10
            δn  = q(k -> f0(k) * (k^2 * Σ / 3 + c.d))
            δΠ  = q(k -> f0(k) / E(k) * (k^4 * Σ / 9 + c.d * k^2 / 3))
            δπr = q(k -> f0(k) / E(k) * k^4 * (4c.Ar - 2(c.Aφ + c.Aη)) / 45)
            δπφ = q(k -> f0(k) / E(k) * k^4 * (4c.Aφ - 2(c.Ar + c.Aη)) / 45)
            δπη = q(k -> f0(k) / E(k) * k^4 * (4c.Aη - 2(c.Ar + c.Aφ)) / 45)
            @test abs(δn) < 1e-10 * n
            @test δΠ ≈ piB_Q rtol = 1e-9
            @test δπφ ≈ piphiphi_Q rtol = 1e-9
            @test δπη ≈ pietaeta_Q rtol = 1e-9
            @test δπr ≈ -(piphiphi_Q + pietaeta_Q) rtol = 1e-9
        end
    end

    @testset "charm through a surface element: δf₁ adds the flux of ν, δf₂ adds nothing" begin
        rng = MersenneTwister(7)
        for _ in 1:6
            s = df_node(rng); ut = sqrt(1 + s.ur^2)
            N(f) = quadgk(pt -> 2π * pt * f(pt), 0, 25; rtol = 1e-12)[1]
            N0 = N(pt -> df_f0(s, pt))
            # dΣ_μ V^μ ∝ dτ V^r − dr V^τ; ν^τ = u^r ν^r/u^τ
            flux = (s.dta * s.νr - s.dra * s.ur * s.νr / ut) / (s.n * (s.dta * s.ur - s.dra * ut))
            @test N(pt -> df_f0(s, pt; ν = s.νr / s.n)) / N0 - 1 ≈ flux rtol = 1e-8
            @test abs(N(pt -> df_d2(s, pt))) < 1e-9 * abs(N0)
        end
    end

    # a freeze-out surface in the two_way_1d layout; nur, piphiphi_Q, pietaeta_Q, piB_Q scaled by c
    x_fn((s,)) = SVector(1 + 8s, 7 * (1 - s^2))
    fields(c) = ((s,),) -> begin
        T = 0.156; ur = 0.9 * s * (2 - s); α = -2.0
        n = F.charm_eos(T, α, DF_CH).n
        SVector(T, ur, 0.0, 0.0, 0.0, α, c * 0.2n * (1 + s), -c * 0.04n * T, c * (0.04 - 0.12 * (1 - s)) * n * T, c * 0.06n * T * s)
    end
    fo(c) = FreezeOutResult(F.SplineInterp(F.PieceWiseFunction(x_fn, 0.0, 1.0), (200,)),
                            F.SplineInterp(F.PieceWiseFunction(fields(c), 0.0, 1.0), (200,)))

    @testset "the integrated spectrum carries the same δf₁" begin
        # spectra(pT, fo, part, eos, transport; delta_f) integrates f over η_s and ϕ_s numerically
        x, phi = fo(1.0); eos = F.Heavy_Quark(hadron_list = DF_CH); part = DF_CH.particle_list[1]
        on, off = CharmDeltaF(DF_CH; second_moment = false), nothing
        for a in (0.2, 0.5, 0.8), pt in (0.5, 2.0)
            V(d) = quadgk(ηs -> quadgk(ϕs -> F._dn_pTdpTdx_internal_spectra(pt, 0.0, 0.0, a, ϕs, ηs, x, phi, part, eos,
                                                                           (nothing, nothing); delta_f = d),
                                       0, 2π; rtol = 1e-11)[1], -6, 6; rtol = 1e-11)[1]
            I(d) = F._pointwise_spectra_internal(pt, DF_M, a, x, phi; deg = DF_DEG, delta_f = d)
            @test V(true) / V(false) ≈ I(on) / I(off) rtol = 1e-8
        end
    end

    @testset "spectra_internal on a freeze-out surface" begin
        kw = (deg = DF_DEG, pt_max = 6.0, step = 25)
        plain = F.spectra_internal(DF_M, fo(0.0); kw...)
        # off by default: the fields' nur, piphiphi_Q, pietaeta_Q, piB_Q are ignored
        @test F.spectra_internal(DF_M, fo(1.0); kw...) == plain
        # and exactly zero when they all vanish
        @test F.spectra_internal(DF_M, fo(0.0); kw..., delta_f = CharmDeltaF(DF_CH)) == plain
        # the second moment changes the shape but not dN/dy
        x, phi = fo(1.0)
        lb, rb = F.leftbounds(x)[1], F.rightbounds(x)[1]
        S(pt, d) = quadgk(a -> F._pointwise_spectra_internal(pt, DF_M, a, x, phi; deg = DF_DEG, delta_f = d), lb, rb;
                          rtol = 1e-11)[1]
        d1, d12 = CharmDeltaF(DF_CH; second_moment = false), CharmDeltaF(DF_CH)
        @test abs(S(2.0, d12) / S(2.0, d1) - 1) > 1e-3
        N(d) = quadgk(pt -> 2π * pt * S(pt, d), 0, 20; rtol = 1e-9)[1]
        @test N(d12) ≈ N(d1) rtol = 1e-7
    end
end
