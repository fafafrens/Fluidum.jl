using Fluidum
using Test
using StaticArrays: SVector
using QuadGK: quadgk
using Interpolations: linear_interpolation, Flat
using Bessels: besselk1, besseli0

# A toy conserved charge for `fo_integral`, which reads the density as
# thermodynamic(T, α, eos.hadron_list).pressure. Here n(T, α) = e^α T³.
struct FOToyCharge end
Fluidum.thermodynamic(T, α, ::FOToyCharge) = (pressure = exp(α) * T^3,)
const fo_toy_eos = (hadron_list = FOToyCharge(),)

# Temperature profiles, T(τ, r) in GeV with τ, r in fm.
# A Bjorken-like fireball that is hottest at the centre.
fo_T_mono(τ, r) = 0.6 * (0.4 / τ)^(1 / 3) * exp(-r^2 / (2 * 4.0^2))
# A Gaussian ring whose radius grows after τ = 2: the centre cools through T_fo while a hot
# shell survives, so the isotherm has two radial crossings at the same τ.
fo_T_shell(τ, r) = 0.5 * (0.4 / τ)^(1 / 3) * exp(-(r - (τ < 2.0 ? 0.0 : 0.8 * (τ - 2.0)))^2 / (2 * 1.5^2))

# A conserved current from a stream function ψ(τ, r) = Q (1 - exp(-r²/2w²)), w = w₀ (τ/τ₀)^b:
# τ r N^τ = ∂_r ψ and τ r N^r = -∂_τ ψ satisfy ∂_τ(τ r N^τ) + ∂_r(τ r N^r) = 0 identically.
# The charge on the τ₀ slice up to R is 2π ψ(τ₀, R). Given a flow u^r, the density and the
# diffusion current follow from N^τ = n u^τ + u^r ν/u^τ and N^r = n u^r + ν.
const fo_tau0 = 0.4
fo_w(τ) = 2.0 * (τ / fo_tau0)^0.25
fo_psi(τ, r) = 1 - exp(-r^2 / (2 * fo_w(τ)^2))
fo_ur(τ, r) = 0.2 * tanh(r / 3) * (1 - fo_tau0 / τ)
function fo_toy_fields(Tfun)
    return function (τ, r)
        w = fo_w(τ); dw = 0.25 * w / τ
        g = exp(-r^2 / (2 * w^2))
        Nt = g / (τ * w^2)
        Nr = g * r * dw / (τ * w^3)
        ur = fo_ur(τ, r); ut = sqrt(1 + ur^2)
        n = ut * Nt - ur * Nr
        ν = Nr - n * ur
        T = Tfun(τ, r)
        return SVector(T, ur, 0.0, 0.0, 0.0, log(n / T^3), ν)
    end
end

@testset "freeze-out contour" begin

@testset "isotherm with a hot shell" begin
    Tfo = 0.156
    taus = collect(0.4:0.025:14.0)
    o = Fluidum.FOContourOpts(Tfo = Tfo, rmax = 12.0)

    b = Fluidum.isotherm_branches(fo_T_shell, taus, o)
    @test !isempty(b.outer_tau)
    @test length(b.inner_tau) >= o.min_inner_run
    # the shell freezes out when its peak, T_max(τ) = 0.5 (0.4/τ)^(1/3), reaches Tfo
    @test b.tau_pinch < last(taus)
    @test b.tau_pinch ≈ 0.4 * (0.5 / Tfo)^3 atol = 0.1
    @test maximum(abs.(fo_T_shell.(b.outer_tau, b.outer_r) .- Tfo)) < 1e-3
    @test maximum(abs.(fo_T_shell.(b.inner_tau, b.inner_r) .- Tfo)) < 1e-3

    c = Fluidum.isotherm_contour(fo_T_shell, taus, o)
    @test length(c.lambda) == length(c.tau) == length(c.r) == length(c.branch)
    @test all(diff(c.lambda) .>= 0) && c.lambda[end] > 0
    @test count(==(:outer), c.branch) > 0 && count(==(:inner), c.branch) > 0
    # the outer branch rises in τ, the inner one comes back down, and the path closes at the axis
    @test issorted(c.tau[c.branch .== :outer])
    @test issorted(c.tau[c.branch .== :inner]; rev = true)
    @test c.r[end] == 0.0
    @test c.tau[end] ≈ Fluidum.axis_freeze_time(fo_T_shell, taus, o)

    # keyword form and the open path
    c_open = Fluidum.isotherm_contour(fo_T_shell, taus; Tfo = Tfo, rmax = 12.0, close_axis = false)
    @test length(c_open.tau) == length(c.tau) - 1 && c_open.r[end] > 0

    # an inner branch shorter than min_inner_run is discarded as discretisation noise
    c_short = @test_logs (:warn, r"discarding") Fluidum.isotherm_contour(fo_T_shell, taus;
        Tfo = Tfo, rmax = 12.0, min_inner_run = 10_000)
    @test count(==(:inner), c_short.branch) == 0
    @test c_short.lambda[end] < c.lambda[end]

    # a grid that ends while the shell is still hot leaves the surface open, and says so
    @test_logs (:warn, r"open") Fluidum.isotherm_branches(fo_T_shell, collect(0.4:0.025:12.0), o)
end

@testset "monotone fireball: one outer branch" begin
    Tfo = 0.156
    taus = collect(0.4:0.025:30.0)
    o = Fluidum.FOContourOpts(Tfo = Tfo, rmax = 12.0)
    b = Fluidum.isotherm_branches(fo_T_mono, taus, o)
    @test isempty(b.inner_tau)
    # r_out(τ) against the closed form 4 √(2 ln(T_axis(τ)/Tfo)), to the radial scan resolution
    r_exact = [4.0 * sqrt(2 * log(fo_T_mono(τ, 0.0) / Tfo)) for τ in b.outer_tau]
    @test maximum(abs.(b.outer_r .- r_exact)) < 2 * o.rmax / o.nr

    c = Fluidum.isotherm_contour(fo_T_mono, taus, o)
    @test all(==(:outer), c.branch)
    @test issorted(c.tau)
    # closed at the axis when T(τ, 0) = Tfo, i.e. τ = τ₀ (0.6/Tfo)³
    @test c.r[end] == 0.0
    @test c.tau[end] ≈ 0.4 * (0.6 / Tfo)^3 rtol = 1e-4
end

@testset "τ₀ cap against its closed form" begin
    τ0, r_edge, r_max = 0.4, 4.0, 12.0
    Tc(r) = 0.15 * exp(-(r - r_edge) / 6)
    urc(r) = 0.3 * r / (r + 2)
    αc(r) = -1.0 - 0.1 * r
    fields(τ, r) = SVector(Tc(r), urc(r), 0.0, 0.0, 0.0, αc(r), 0.0)
    cap = Fluidum.freezeout_cap(τ0, r_edge, r_max, fields)
    x, phi = cap

    # a purely τ-like normal, running inwards
    for s in (0.0, 0.3, 0.77, 1.0)
        @test x(s) ≈ SVector(τ0, r_max - s * (r_max - r_edge))
        dτ, dr = Fluidum.jacobian(x, s)
        @test abs(dτ) < 1e-10
        @test dr ≈ -(r_max - r_edge)
    end
    @test_throws ArgumentError Fluidum.freezeout_cap(τ0, r_max, r_edge, fields)

    m, deg = 1.87, 2
    pts = range(0.5, 4.0, 8)
    fmGeV = 5.068   # the constant used inside internal_thermal_spectra
    closed(pt) = quadgk(r_edge, r_max; rtol = 1e-10) do r
        T, ur, α = Tc(r), urc(r), αc(r)
        mt = sqrt(m^2 + pt^2)
        (τ0 * r / (2π^2)) * mt * besselk1(mt * sqrt(1 + ur^2) / T) * besseli0(pt * ur / T) *
            exp(α) * fmGeV^3 * deg
    end[1]
    S = getindex.(Fluidum.spectra_internal(m, cap; pt_min = 0.5, pt_max = 4.0, step = 8, deg = deg), 1)
    @test all(S .> 0)
    @test S ≈ closed.(pts) rtol = 1e-5
end

@testset "arc length against a τ-parameterised surface" begin
    # The outer edge of the monotone fireball up to τ₁, before it turns towards the axis,
    # parameterised by τ as a single-valued r_fo(τ) and by arc length.
    Tfo = 0.156
    rfo(τ) = 4.0 * sqrt(2 * log(fo_T_mono(τ, 0.0) / Tfo))
    τ0, τ1 = 0.4, 10.0
    fields = fo_toy_fields(fo_T_mono)
    by_tau = FreezeOutResult(
        SplineInterp(Fluidum.PieceWiseFunction(t -> SVector(t[1], rfo(t[1])), τ0, τ1), (300,)),
        SplineInterp(Fluidum.PieceWiseFunction(t -> fields(t[1], rfo(t[1])), τ0, τ1), (300,)))
    τs = collect(range(τ0, τ1, 400))
    by_arc = Fluidum.freezeout_from_path(τs, rfo.(τs), fields)

    I_tau = fo_integral(by_tau, fo_toy_eos)[1]
    I_arc = fo_integral(by_arc, fo_toy_eos)[1]
    @test I_arc ≈ I_tau rtol = 1e-4

    S_tau = getindex.(Fluidum.spectra_internal(1.87, by_tau; pt_min = 0.5, pt_max = 4.0, step = 6), 1)
    S_arc = getindex.(Fluidum.spectra_internal(1.87, by_arc; pt_min = 0.5, pt_max = 4.0, step = 6), 1)
    @test S_arc ≈ S_tau rtol = 1e-4

    # sampled fields give the same surface as the callable
    by_samples = Fluidum.freezeout_from_path(τs, rfo.(τs), fields.(τs, rfo.(τs)))
    @test fo_integral(by_samples, fo_toy_eos)[1] ≈ I_arc rtol = 1e-4

    # the orientation is the order of the points: reversing the path flips the flux
    rev = Fluidum.freezeout_from_path(reverse(τs), reverse(rfo.(τs)), fields)
    @test fo_integral(rev, fo_toy_eos)[1] ≈ -I_arc rtol = 1e-6
end

@testset "charge conservation: contour + cap = charge at τ₀" begin
    Tfo, R = 0.156, 12.0
    for (Tfun, taus) in ((fo_T_mono, collect(0.4:0.025:30.0)), (fo_T_shell, collect(0.4:0.025:14.0)))
        fields = fo_toy_fields(Tfun)
        fo = Fluidum.freezeout_from_isotherm(Tfun, fields, taus; Tfo = Tfo, rmax = R, npts = 600)
        @test fo.cap !== nothing
        I_surface = fo_integral(fo.surface, fo_toy_eos; rtol = 1e-9)[1]
        I_cap = fo_integral(fo.cap, fo_toy_eos; rtol = 1e-9)[1]
        Q0 = 2π * fo_psi(fo_tau0, R)
        @test I_surface > 0 && I_cap > 0
        @test I_surface + I_cap ≈ Q0 rtol = 1e-4
    end
    # without the segment that closes it at the axis, the monotone surface misses charge
    fields = fo_toy_fields(fo_T_mono)
    fo = Fluidum.freezeout_from_isotherm(fo_T_mono, fields, collect(0.4:0.025:30.0);
        Tfo = Tfo, rmax = R, npts = 600, close_axis = false)
    I_open = fo_integral(fo.surface, fo_toy_eos; rtol = 1e-9)[1] + fo_integral(fo.cap, fo_toy_eos; rtol = 1e-9)[1]
    @test !isapprox(I_open, 2π * fo_psi(fo_tau0, R); rtol = 1e-4)
end

@testset "Fluidum hydro: freeze_out_routine against the isotherm contour" begin
    eos = Heavy_Quark()
    params = Fluidum.FluidProperties(eos, ZeroViscosity(), ZeroBulkViscosity(), ZeroDiffusion())
    disc = CartesianDiscretization(OriginInterval(100, 12.0))
    df = DiscreteFields(Fluidum.HQ_viscous_1d(), disc, Float64)
    Tini(r) = 0.4 / (exp((r - 4) / 0.8) + 1) + 0.05
    phi = set_array(Tini, :temperature, df)
    set_array!(phi, r -> -3.0, :α, df)
    tspan = (0.4, 15.0)
    Tfo = 0.156

    fo_up = freeze_out_routine(df, Fluidum.matrix1d_visc_HQ!, params, phi, tspan; Tfo = Tfo)
    I_up = fo_integral(fo_up, eos)[1]
    S_up = getindex.(Fluidum.spectra_internal(1.87, fo_up; pt_min = 0.5, pt_max = 4.0, step = 6), 1)

    # 1. the same surface, resampled and parameterised by arc length
    x, f = fo_up
    as = range(Fluidum.leftbounds(x)[1], Fluidum.rightbounds(x)[1], 400)
    pts = x.(as)
    fo_rs = Fluidum.freezeout_from_path(getindex.(pts, 1), getindex.(pts, 2), f.(as))
    @test fo_integral(fo_rs, eos)[1] ≈ I_up rtol = 1e-4
    S_rs = getindex.(Fluidum.spectra_internal(1.87, fo_rs; pt_min = 0.5, pt_max = 4.0, step = 6), 1)
    @test S_rs ≈ S_up rtol = 1e-4

    # 2. the isotherm of the same hydro solution, found by isotherm_contour
    sol = Fluidum.oneshoot(df, Fluidum.matrix1d_visc_HQ!, params, phi, tspan)
    rg = [disc.grid[i][1] for i in 2:length(disc.grid)-1]
    taus = collect(0.4:0.02:15.0)
    tab = [SVector{7}(sol(τ)[:, i+1]) for τ in taus, i in eachindex(rg)]
    itp = linear_interpolation((taus, rg), tab; extrapolation_bc = Flat())
    fields(τ, r) = itp(τ, r)
    T_hydro(τ, r) = itp(τ, r)[1]
    fo = Fluidum.freezeout_from_isotherm(T_hydro, fields, taus; Tfo = Tfo, rmax = rg[end], nr = 2000)
    @test all(==(:outer), fo.contour.branch)
    @test first(fo.contour.r) ≈ first(x(Fluidum.leftbounds(x)[1])[2]) rtol = 1e-3
    @test fo_integral(fo.surface, eos)[1] ≈ I_up rtol = 2e-3
    S_new = getindex.(Fluidum.spectra_internal(1.87, fo.surface; pt_min = 0.5, pt_max = 4.0, step = 6), 1)
    @test S_new ≈ S_up rtol = 2e-3
    # the cold exterior at τ₀ adds a small, positive amount
    @test 0 < fo_integral(fo.cap, eos)[1] < 1e-3 * I_up
end

end
