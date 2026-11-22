# test/test_2p1d_viscous.jl -- the 2+1D viscous solver, `matrix2d_visc!`.
#
# Checked here: closed-form characteristic speeds, isotropy of the shear rows, and the Bjorken limit.
using Test, LinearAlgebra
# no Printf: it is not in the `Pkg.test` environment.
using Fluidum, OrdinaryDiffEq
const _F = Fluidum

_prs(T, eos)  = Fluidum.pressure(T, eos)
_ent(T, eos)  = Fluidum.pressure_derivative(T, Val(1), eos)
_dst(T, eos)  = Fluidum.pressure_derivative(T, Val(2), eos)
_ene(T, eos)  = T*_ent(T, eos) - _prs(T, eos)
_enth(T, eos) = T*_ent(T, eos)
_cs2(T, eos)  = _ent(T, eos)/(T*_dst(T, eos))

const _EOS  = FluiduMEoS()
const _PAR  = Fluidum.FluidProperties(_EOS, QGPViscosity(0.1, 0.2),
                                      SimpleBulkViscosity(0.083, 15.0), ZeroDiffusion())
const _PARI = Fluidum.FluidProperties(_EOS, ZeroViscosity(), ZeroBulkViscosity(), ZeroDiffusion())

function _coeffs(T, par = _PAR)
    th = thermodynamic(T, par.eos)
    (Fluidum.viscosity(T, th, par.shear),      Fluidum.τ_shear(T, th, par.shear),
     Fluidum.bulk_viscosity(T, th, par.bulk),  Fluidum.τ_bulk(T, th, par.bulk))
end
function _sys(T, ux, uy, piyy, pizz, pixy, piB, tau, par = _PAR; m = Fluidum.matrix2d_visc!)
    A = [zeros(7,7), zeros(7,7)]; S = zeros(7)
    m(A, S, [T,ux,uy,piyy,pizz,pixy,piB], tau, (1.0,1.0), par)
    A[1], A[2], S
end
_spec(A) = sort(real.(eigvals(A)))

@testset "2+1D viscous -- field layout" begin
    f = Fluidum.viscous_2d()
    @test length(f.field) == 7
    @test collect(Fluidum.names(f)) == [:temperature, :ux, :uy, :piyy, :pizz, :pixy, :piB]
    @test Fluidum.get_index(:pizz, f) == 5      # the mixed pi^eta_eta -- see the Bjorken testset
end

@testset "2+1D viscous -- shear closure" begin
    # the three stored components must close onto a traceless, u-orthogonal pi^{mu nu} through
    #   pi^{xx} = [ -(pi^{yy}+pi^{zz})(u^t)^2 + 2 u^x u^y pi^{xy} + (u^y)^2 pi^{yy} ] / (1 + (u^y)^2)
    g = Diagonal([-1.0, 1.0, 1.0, 1.0])
    for ux in (-2.0, 0.0, 0.7), uy in (-1.5, 0.0, 2.0),
        (piyy, pizz, pixy) in ((0.1,-0.2,0.05), (-0.03,0.07,-0.11))
        ut  = sqrt(1 + ux^2 + uy^2)
        pixx = (-(piyy + pizz)*ut^2 + 2ux*uy*pixy + uy^2*piyy)/(1 + uy^2)
        pitx = (ux*pixx + uy*pixy)/ut; pity = (ux*pixy + uy*piyy)/ut
        pitt = (ux*pitx + uy*pity)/ut
        P = [pitt pitx pity 0.0; pitx pixx pixy 0.0; pity pixy piyy 0.0; 0.0 0.0 0.0 pizz]
        u = [ut, ux, uy, 0.0]; sc = maximum(abs, P)
        @test abs(tr(g*P))/sc < 1e-13
        @test maximum(abs, (u'*g*P)')/sc < 1e-13
    end
end

@testset "2+1D viscous -- characteristic speeds against closed forms" begin
    for T in (0.15, 0.30, 0.60)
        eta, taupi, zeta, tauPi = _coeffs(T); w = _enth(T, _EOS)
        # the Israel-Stewart front, at rest
        vf = sqrt(_cs2(T,_EOS) + (4/3)*eta/(w*taupi) + zeta/(w*tauPi))
        @test maximum(abs, _spec(_sys(T,0,0,0,0,0,0,2.0)[1])) ≈ vf rtol=1e-12
        # The transverse shear channel. Israel-Stewart gives sqrt(eta/(w tau_pi)) = sqrt(C_s),
        # independent of T and the EoS. The legacy matrix (one-sided sigma^{xy}) gives sqrt(2 C_s)
        # along x and no channel along y; that value is checked below with the legacy matrix.
        @test minimum(abs.(_spec(_sys(T,0,0,0,0,0,0,2.0)[1]) .- sqrt(0.2))) < 1e-12
        @test minimum(abs.(_spec(_sys(T,0,0,0,0,0,0,2.0; m = Fluidum.matrix2d_visc!(; legacy = true))[1]) .-
                           sqrt(0.4))) < 1e-12
        # with viscosity off, the extreme characteristics are the relativistic sound cone
        for U in (0.0, 1.2)
            c = sqrt(_cs2(T,_EOS)); v = U/sqrt(1+U^2)
            ev = _spec(_sys(T,U,0,0,0,0,0,2.0,_PARI)[1])
            @test ev[end] ≈ (v+c)/(1+v*c) rtol=1e-11
            @test ev[1]   ≈ (v-c)/(1-v*c) rtol=1e-11
        end
        # A_t must be invertible everywhere the solver goes
        @test isfinite(sum(_sys(T,0.3,-0.2,0.01,-0.02,0.003,0.005,2.0)[1]))
    end
end

# Isotropy and the symmetric shear gradient. The derived matrix is the default; the legacy matrix
# (`matrix2d_visc!(; legacy = true)`) fails both checks.
@testset "2+1D viscous -- isotropy and the symmetric sigma^{xy} (derived matrix)" begin
    for T in (0.20, 0.45)
        Ax, Ay, _ = _sys(T, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0)
        # at rest and isotropic, the x- and y-flux spectra must coincide
        @test maximum(abs, _spec(Ax) .- _spec(Ay)) < 1e-11
        # and pi^{xy} must be driven by the symmetric gradient: eta/tau_pi on each half
        eta, taupi, _, _ = _coeffs(T)
        @test Ax[6,3] ≈ eta/taupi rtol=1e-10
        @test Ay[6,2] ≈ eta/taupi rtol=1e-10
        # the legacy matrix: anisotropic at rest
        Lx, Ly, _ = _sys(T, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0; m = Fluidum.matrix2d_visc!(; legacy = true))
        @test maximum(abs, _spec(Lx) .- _spec(Ly)) > 1e-3
    end
    # the closure at its default is the positional method
    ϕ = [0.3, 0.2, -0.1, 0.01, -0.02, 0.003, 0.004]
    A = [zeros(7,7), zeros(7,7)]; S = zeros(7); A2 = [zeros(7,7), zeros(7,7)]; S2 = zeros(7)
    Fluidum.matrix2d_visc!(A, S, ϕ, 1.1, (1.0, 1.0), _PAR)
    Fluidum.matrix2d_visc!()(A2, S2, ϕ, 1.1, (1.0, 1.0), _PAR)
    @test A == A2 && S == S2
end

@testset "2+1D viscous -- Bjorken is the 0+1D Israel-Stewart (+ delta_pipi) system" begin
    # A transversely uniform state has no transverse gradient, so the 2-D solver must integrate
    #   De = -(e + P + Pi + pi^eta_eta)/tau,
    #   tau_pi D pi^eta_eta = -pi^eta_eta - 4 eta/(3 tau) - delta_pipi pi^eta_eta/tau,
    #   tau_Pi D Pi = -Pi - zeta/tau
    # written here from the physics rather than from the matrix (theta = 1/tau on Bjorken).
    # This pins the convention: pi^{zz} is the
    # mixed pi^eta_eta = tau^2 pi^{eta eta}. delta_pipi = (4/3) tau_pi by default and 0 for
    # `matrix2d_visc!(; delta_pipi = false)`. Each value is checked against its own ODE and must miss
    # the other's by >5 %.
    tau0, tau1, T0 = 0.4, 3.0, 0.35
    function bjref(fac)
        function rhs!(du, U, p, tau)
            T, pe, PB = U
            eta, taupi, zeta, tauPi = _coeffs(T)
            du[1] = (-(_ene(T,_EOS) + _prs(T,_EOS) + PB + pe)/tau)/(T*_dst(T,_EOS))
            du[2] = (-pe - 4eta/(3tau) - fac*(4/3)*taupi*pe/tau)/taupi
            du[3] = (-PB - zeta/tau)/tauPi
        end
        solve(ODEProblem(rhs!, [T0,0.0,0.0], (tau0,tau1)), Tsit5(); reltol=1e-12, abstol=1e-14)(tau1)
    end

    d = DiscreteFields(Fluidum.viscous_2d(),
        CartesianDiscretization(Fluidum.SymmetricInterval(24, 12.0),
                                Fluidum.SymmetricInterval(24, 12.0)), Float64)
    ic = 13
    for fac in (1.0, 0.0)
        # the default positional method for 4/3 tau_pi, the keyword closure for plain IS
        m = fac == 1.0 ? Fluidum.matrix2d_visc! : Fluidum.matrix2d_visc!(; delta_pipi = false)
        ref = bjref(fac)
        phi = set_array((x,y)->T0, :temperature, d)
        for fl in (:ux,:uy,:piyy,:pizz,:pixy,:piB); set_array!(phi, (x,y)->0.0, fl, d); end
        u = Fluidum.oneshoot(d, m, _PAR, phi, (tau0,tau1);
                             reltol=1e-11, abstol=1e-13)(tau1)
        @test u[1,ic,ic] ≈ ref[1] rtol=1e-10
        @test u[5,ic,ic] ≈ ref[2] rtol=1e-10
        @test u[7,ic,ic] ≈ ref[3] rtol=1e-10
        @test abs(u[5,ic,ic]/bjref(1.0 - fac)[2] - 1) > 0.05
        # tracelessness is dynamical: pi^{yy} and pi^{zz} are advanced by separate rows
        @test u[4,ic,ic] ≈ -u[5,ic,ic]/2 rtol=1e-11
        # a uniform state must stay uniform and at rest
        @test maximum(abs, u[1,2:25,2:25] .- u[1,ic,ic])/u[1,ic,ic] < 1e-13
        @test maximum(abs, u[2,2:25,2:25]) < 1e-13
        @test maximum(abs, u[3,2:25,2:25]) < 1e-13
    end

    # and with viscosity off, tau*s is conserved on the Bjorken solution
    phi = set_array((x,y)->T0, :temperature, d)
    for f in (:ux,:uy,:piyy,:pizz,:pixy,:piB); set_array!(phi, (x,y)->0.0, f, d); end
    ui = Fluidum.oneshoot(d, Fluidum.matrix2d_visc!, _PARI, phi, (tau0,tau1);
                          reltol=1e-11, abstol=1e-13)(tau1)
    @test tau1*_ent(ui[1,ic,ic], _EOS) ≈ tau0*_ent(T0, _EOS) rtol=1e-9
end

@testset "2+1D viscous -- the 1-D and 2-D transport paths agree" begin
    # `matrix2d_visc!` passes a whole `Thermodynamic` and `matrix1d_visc!` a scalar entropy.
    # The formulas are identical and must give the same numbers.
    for T in (0.15, 0.25, 0.35, 0.5)
        th = thermodynamic(T, _EOS); dpt = _ent(T,_EOS); dptt = _dst(T,_EOS)
        @test Fluidum.viscosity(T, dpt, _PAR.shear)     ≈ Fluidum.viscosity(T, th, _PAR.shear)     rtol=1e-12
        @test Fluidum.bulk_viscosity(T, dpt, _PAR.bulk) ≈ Fluidum.bulk_viscosity(T, th, _PAR.bulk) rtol=1e-12
    end
end
