# the convention here are α, nux, nuy on a prescribed background
# background_fields(tau, x, y) = (T, ux, uy, dtT, dxT, dyT, dtux, dxux, dyux, dtuy, dxuy, dyuy)
# Fick closure:       Nabla.(n u + nu) = 0, tau_n Delta D nu + nu + kappa Delta grad(alpha) = 0
# consistent closure: Nabla.(n u + nu) = 0, (D_s/T) Delta^mu_l Nabla_m T_Q^{ml} + nu^mu = 0, P0 = nT
# matrix2d_visc_HQ_BG!: T, ux, uy, piyy, pizz, pixy, piB, α, nux, nuy, 2+1D viscous hydro with passive charm (Fick)

"""
    matrix2d_HQ_BG!(; background_fields, dmn_eps = 1e-6, α_max = 200.0)
Charm first moment (α, ν^x, ν^y) on a prescribed 2+1D background, Fick closure.
"""
function matrix2d_HQ_BG!(; dmn_eps = 1e-6, background_fields, α_max = 200.0)
    (A_i, Source, ϕ, tau, X, params) -> matrix2d_HQ_BG!(A_i, Source, ϕ, tau, X, params;
        dmn_eps = dmn_eps, background_fields = background_fields, α_max = α_max)
end

function matrix2d_HQ_BG!(A_i, Source, ϕ, tau, X, params;
        dmn_eps = 1e-6, background_fields = nothing, α_max = 200.0)

    T, ux, uy, dtT, dxT, dyT, dtux, dxux, dyux, dtuy, dxuy, dyuy = background_fields(tau, X[1], X[2])

    α_safe = clamp(ϕ[1], -α_max, α_max)
    thermo = thermodynamic(T, α_safe, params.eos.hadron_list)
    n = thermo.pressure
    dn_dT, dn_da = thermo.pressure_derivative
    dn_da += dmn_eps        # keeps inv(At) regular where n -> 0 (the α row's A_t entry is u^tau dn_da)
    κ    = diffusion_hadron(T, α_safe, params.eos, params.diffusion)
    taun = τ_diffusion_hadron(T, α_safe, params.eos, params.diffusion)

    At, Ax, Ay, src = two_d_HQ_BG_matrix(α_safe, ϕ[2], ϕ[3], tau, T, ux, uy,
                                         dtT, dxT, dyT, dtux, dxux, dyux, dtuy, dxuy, dyuy,
                                         n, dn_dT, dn_da, taun, κ)

    Ainv = inv(At)
    A_mul_B!(A_i[1], Ainv, Ax)
    A_mul_B!(A_i[2], Ainv, Ay)
    jgemvavx!(Source, Ainv, src)
end

"""
    matrix2d_HQ_BG_consistent!(; background_fields, dmn_eps = 1e-6, α_max = 200.0)
Charm first moment (α, ν^x, ν^y) on a prescribed 2+1D background, consistent closure:
at rest the current relaxes to -κ∇α - (D_s/T)(n + T ∂_T n)∇T.
"""
function matrix2d_HQ_BG_consistent!(; dmn_eps = 1e-6, background_fields, α_max = 200.0)
    (A_i, Source, ϕ, tau, X, params) -> matrix2d_HQ_BG_consistent!(A_i, Source, ϕ, tau, X, params;
        dmn_eps = dmn_eps, background_fields = background_fields, α_max = α_max)
end

function matrix2d_HQ_BG_consistent!(A_i, Source, ϕ, tau, X, params;
        dmn_eps = 1e-6, background_fields = nothing, α_max = 200.0)

    T, ux, uy, dtT, dxT, dyT, dtux, dxux, dyux, dtuy, dxuy, dyuy = background_fields(tau, X[1], X[2])

    ch = charm_eos(T, clamp(ϕ[1], -α_max, α_max), params.eos.hadron_list)
    Ds = DsT(params.diffusion, T) / T / fmGeV

    # dn/dα + dmn_eps keeps inv(At) regular where n -> 0; it enters kappa = D_s dn/dα too
    At, Ax, Ay, src = two_d_HQ_BG_consistent_matrix(ϕ, tau, T, ux, uy, dtT, dxT, dyT, dtux, dxux, dyux, dtuy, dxuy, dyuy,
                                                    ch.n, ch.dndT, ch.dnda + dmn_eps, ch.h, ch.dhdT, Ds)

    Ainv = inv(At)
    A_mul_B!(A_i[1], Ainv, Ax)
    A_mul_B!(A_i[2], Ainv, Ay)
    jgemvavx!(Source, Ainv, src)
end

"""
    matrix2d_visc_HQ_BG!(; dmn_eps = 1e-6, α_max = 200.0, legacy = false, delta_pipi = true)
2+1D viscous hydro with passive charm (Fick closure), 10 fields (`HQ_viscous_2d`); `legacy` and
`delta_pipi` as in `matrix2d_visc!`.
"""
function matrix2d_visc_HQ_BG!(; dmn_eps = 1e-6, α_max = 200.0, legacy::Bool = false, delta_pipi::Bool = true)
    (A_i, Source, ϕ, tau, X, params) -> matrix2d_visc_HQ_BG!(A_i, Source, ϕ, tau, X, params;
        dmn_eps = dmn_eps, α_max = α_max, legacy = legacy, delta_pipi = delta_pipi)
end

function matrix2d_visc_HQ_BG!(A_i, Source, ϕ, tau, X, params; dmn_eps = 1e-6, α_max = 200.0,
                              legacy::Bool = false, delta_pipi::Bool = true)
    T  = ϕ[1]
    p    = pressure(T, params.eos)
    dpt  = pressure_derivative(T, Val(1), params.eos)     # entropy
    dptt = pressure_derivative(T, Val(2), params.eos)
    η  = viscosity(T, dpt, params.shear)
    τS = τ_shear(T, dpt, params.shear)
    ζ  = bulk_viscosity(T, dpt, params.bulk)
    τB = τ_bulk(T, dpt, dptt, params.bulk)

    # rows 1-7: the medium, as in matrix2d_visc!
    uh = SVector{7}(ϕ[1], ϕ[2], ϕ[3], ϕ[4], ϕ[5], ϕ[6], ϕ[7])
    Ath, Axh, Ayh, srch = legacy ? two_d_viscous_matrix(uh, tau, p, dpt, dptt, ζ, η, τS, τB) :
                       delta_pipi ? two_d_viscous_matrix_derived(uh, tau, p, dpt, dptt, ζ, η, τS, τB) :
                                    two_d_viscous_matrix_derived(uh, tau, p, dpt, dptt, ζ, η, τS, τB, zero(τS))

    # rows 8-10: the charm, in the columns (T, ux, uy, α, nux, nuy) = (1, 2, 3, 8, 9, 10)
    α_safe = clamp(ϕ[8], -α_max, α_max)
    thermo = thermodynamic(T, α_safe, params.eos.hadron_list)
    n = thermo.pressure
    dn_dT, dn_da = thermo.pressure_derivative
    dn_da += dmn_eps        # keeps inv(At) regular where n -> 0 (the α row's A_t entry is u^tau dn_da)
    κ    = diffusion_hadron(T, α_safe, params.eos, params.diffusion)
    taun = τ_diffusion_hadron(T, α_safe, params.eos, params.diffusion)
    AtL, AxL, AyL, srcL = two_d_HQ_live_rows(α_safe, ϕ[9], ϕ[10], tau, T, ϕ[2], ϕ[3],
                                             n, dn_dT, dn_da, taun, κ)

    At = zeros(MMatrix{10,10,Float64}); Ax = zeros(MMatrix{10,10,Float64})
    Ay = zeros(MMatrix{10,10,Float64}); src = zeros(MVector{10,Float64})
    @inbounds for j in 1:7, i in 1:7          # passive: (1:7, 8:10) stays zero
        At[i,j] = Ath[i,j]; Ax[i,j] = Axh[i,j]; Ay[i,j] = Ayh[i,j]
    end
    @inbounds for i in 1:7; src[i] = srch[i]; end
    cols = (1, 2, 3, 8, 9, 10)
    @inbounds for k in 1:6, i in 1:3
        At[7+i, cols[k]] = AtL[i, k]; Ax[7+i, cols[k]] = AxL[i, k]; Ay[7+i, cols[k]] = AyL[i, k]
    end
    @inbounds for i in 1:3; src[7+i] = srcL[i]; end

    Ainv = inv(SMatrix{10,10}(At))
    A_mul_B!(A_i[1], Ainv, SMatrix{10,10}(Ax))
    A_mul_B!(A_i[2], Ainv, SMatrix{10,10}(Ay))
    jgemvavx!(Source, Ainv, SVector{10}(src))
end


# generated by symbolic computation (Mathematica)
"""
    two_d_HQ_BG_matrix(α, nux, nuy, tau, T, ux, uy, dtT, dxT, dyT,
                       dtux, dxux, dyux, dtuy, dxuy, dyuy, n, dn_dT, dn_da, taun, kappa)
(A_t, A_x, A_y, src) of the charm first moment phi = (alpha, nu^x, nu^y) on a prescribed 2+1D background,
source on the LHS: Nabla_mu (n u^mu + nu^mu) = 0, tau_n Delta^mu_b D nu^b + nu^mu + kappa Delta^{mu b} d_b alpha = 0.
"""
function two_d_HQ_BG_matrix(α, nux, nuy, tau, T, ux, uy, dtT, dxT, dyT,
                            dtux, dxux, dyux, dtuy, dxuy, dyuy,
                            n, dn_dT, dn_da, taun, kappa)
    At = SMatrix{3,3}((
        sqrt(1 + ux^2 + uy^2)*(dn_da),
        kappa*ux*sqrt(1 + ux^2 + uy^2),
        kappa*uy*sqrt(1 + ux^2 + uy^2),
        ux/sqrt(1 + ux^2 + uy^2),
        taun*sqrt(1 + ux^2 + uy^2),
        0,
        uy/sqrt(1 + ux^2 + uy^2),
        0,
        taun*sqrt(1 + ux^2 + uy^2)))

    Ax = SMatrix{3,3}((
        ux*(dn_da),
        kappa*(1 + ux^2),
        kappa*ux*uy,
        1,
        taun*ux,
        0,
        0,
        0,
        taun*ux))

    Ay = SMatrix{3,3}((
        uy*(dn_da),
        kappa*ux*uy,
        kappa*(1 + uy^2),
        0,
        taun*uy,
        0,
        1,
        0,
        taun*uy))

    src = SVector{3}((
        (-(dtux*nux*tau*ux^2) - dtuy*nux*tau*ux*uy - dtux*nuy*tau*ux*uy - dtuy*nuy*tau*uy^2 + dtux*nux*tau*(1 + ux^2 + uy^2) + dtuy*nuy*tau*(1 + ux^2 + uy^2) + nux*ux*(1 + ux^2 + uy^2) + dtux*n*tau*ux*(1 + ux^2 + uy^2) + nuy*uy*(1 + ux^2 + uy^2) + dtuy*n*tau*uy*(1 + ux^2 + uy^2) + dxux*n*tau*(1 + ux^2 + uy^2)^(3/2) + dyuy*n*tau*(1 + ux^2 + uy^2)^(3/2) + n*(1 + ux^2 + uy^2)^2 + dxT*tau*ux*(1 + ux^2 + uy^2)^(3/2)*(dn_dT) + dyT*tau*uy*(1 + ux^2 + uy^2)^(3/2)*(dn_dT) + dtT*tau*(1 + ux^2 + uy^2)^2*(dn_dT))/(tau*(1 + ux^2 + uy^2)^(3/2)),
        (nux*(1 + dxuy*taun*ux^3*uy + uy^2 - taun*ux*(1 + uy^2)*(dyux*uy + dtux*sqrt(1 + ux^2 + uy^2)) + ux^2*(1 + dyuy*taun*uy^2 - dxux*taun*(1 + uy^2) + dtuy*taun*uy*sqrt(1 + ux^2 + uy^2))) - nuy*taun*ux*(dxuy*ux^3 + dyuy*uy + dtuy*sqrt(1 + ux^2 + uy^2) + ux^2*(-(dxux*uy) + dyuy*uy + dtuy*sqrt(1 + ux^2 + uy^2)) + ux*(dxuy - uy*(dyux*uy + dtux*sqrt(1 + ux^2 + uy^2)))))/(1 + ux^2 + uy^2),
        (nuy*(1 - dxuy*taun*ux^3*uy + (1 - dyuy*taun)*uy^2 - dtuy*taun*uy*sqrt(1 + ux^2 + uy^2) + ux^2*(1 + (dxux - dyuy)*taun*uy^2 - dtuy*taun*uy*sqrt(1 + ux^2 + uy^2)) + taun*ux*uy*(-dxuy + uy*(dyux*uy + dtux*sqrt(1 + ux^2 + uy^2)))) - nux*taun*uy*(dxux*ux + dyux*uy^3 + dtux*sqrt(1 + ux^2 + uy^2) + uy^2*(dxux*ux - dyuy*ux + dtux*sqrt(1 + ux^2 + uy^2)) + uy*(dyux - ux*(dxuy*ux + dtuy*sqrt(1 + ux^2 + uy^2)))))/(1 + ux^2 + uy^2)))

    return (At, Ax, Ay, src)
end


"""
    two_d_HQ_live_rows(α, nux, nuy, tau, T, ux, uy, n, dn_dT, dn_da, taun, kappa)
The same rows with the medium as fields: 3x6 blocks in the columns (T, u^x, u^y, alpha, nu^x, nu^y), and `src`.
"""
function two_d_HQ_live_rows(α, nux, nuy, tau, T, ux, uy, n, dn_dT, dn_da, taun, kappa)
    AtL = SMatrix{3,6}((
        sqrt(1 + ux^2 + uy^2)*(dn_dT),
        0,
        0,
        (-(nuy*ux*uy) + nux*(1 + uy^2) + n*ux*(1 + ux^2 + uy^2))/(1 + ux^2 + uy^2)^(3/2),
        (taun*ux*(nuy*ux*uy - nux*(1 + uy^2)))/sqrt(1 + ux^2 + uy^2),
        -((taun*uy*(nux - nuy*ux*uy + nux*uy^2))/sqrt(1 + ux^2 + uy^2)),
        (nuy*(1 + ux^2) - nux*ux*uy + n*uy*(1 + ux^2 + uy^2))/(1 + ux^2 + uy^2)^(3/2),
        -((taun*ux*(nuy + nuy*ux^2 - nux*ux*uy))/sqrt(1 + ux^2 + uy^2)),
        (taun*uy*(-(nuy*(1 + ux^2)) + nux*ux*uy))/sqrt(1 + ux^2 + uy^2),
        sqrt(1 + ux^2 + uy^2)*(dn_da),
        kappa*ux*sqrt(1 + ux^2 + uy^2),
        kappa*uy*sqrt(1 + ux^2 + uy^2),
        ux/sqrt(1 + ux^2 + uy^2),
        taun*sqrt(1 + ux^2 + uy^2),
        0,
        uy/sqrt(1 + ux^2 + uy^2),
        0,
        taun*sqrt(1 + ux^2 + uy^2)))

    AxL = SMatrix{3,6}((
        ux*(dn_dT),
        0,
        0,
        n,
        (taun*ux^2*(nuy*ux*uy - nux*(1 + uy^2)))/(1 + ux^2 + uy^2),
        (taun*ux*uy*(nuy*ux*uy - nux*(1 + uy^2)))/(1 + ux^2 + uy^2),
        0,
        -((taun*ux^2*(nuy + nuy*ux^2 - nux*ux*uy))/(1 + ux^2 + uy^2)),
        (taun*ux*uy*(-(nuy*(1 + ux^2)) + nux*ux*uy))/(1 + ux^2 + uy^2),
        ux*(dn_da),
        kappa*(1 + ux^2),
        kappa*ux*uy,
        1,
        taun*ux,
        0,
        0,
        0,
        taun*ux))

    AyL = SMatrix{3,6}((
        uy*(dn_dT),
        0,
        0,
        0,
        (taun*ux*uy*(nuy*ux*uy - nux*(1 + uy^2)))/(1 + ux^2 + uy^2),
        -((taun*uy^2*(nux - nuy*ux*uy + nux*uy^2))/(1 + ux^2 + uy^2)),
        n,
        -((taun*ux*uy*(nuy + nuy*ux^2 - nux*ux*uy))/(1 + ux^2 + uy^2)),
        (taun*uy^2*(-(nuy*(1 + ux^2)) + nux*ux*uy))/(1 + ux^2 + uy^2),
        uy*(dn_da),
        kappa*ux*uy,
        kappa*(1 + uy^2),
        0,
        taun*uy,
        0,
        1,
        0,
        taun*uy))

    srcL = SVector{3}((
        (nux*ux + nuy*uy + n*(1 + ux^2 + uy^2))/(tau*sqrt(1 + ux^2 + uy^2)),
        nux,
        nuy))

    return (AtL, AxL, AyL, srcL)
end
# end of generated code


# generated by symbolic computation (Mathematica)
"""
    two_d_HQ_BG_consistent_matrix(u, tau, T, ux, uy, dtT, dxT, dyT, dtux, dxux, dyux, dtuy, dxuy, dyuy, n, dndT, dnda, h, dhdT, Ds)
Charm first moment on a prescribed background, consistent closure; 2+1D Cartesian.
Fields u = (α, nux, nuy); returns (A_t, A_x, A_y, source).
"""
function two_d_HQ_BG_consistent_matrix(u, tau, T, ux, uy, dtT, dxT, dyT, dtux, dxux, dyux, dtuy, dxuy, dyuy, n, dndT, dnda, h, dhdT, Ds)
    α, nux, nuy = u[1], u[2], u[3]
    utau = sqrt(1 + ux^2 + uy^2)
    At = SMatrix{3,3}((
        dnda*utau,
        dnda*Ds*utau*ux,
        dnda*Ds*utau*uy,
        ux/utau,
        (Ds*h*utau)/T,
        0,
        uy/utau,
        0,
        (Ds*h*utau)/T))

    Ax = SMatrix{3,3}((
        dnda*ux,
        dnda*Ds*(1 + ux^2),
        dnda*Ds*ux*uy,
        1,
        (Ds*h*ux)/T,
        0,
        0,
        0,
        (Ds*h*ux)/T))

    Ay = SMatrix{3,3}((
        dnda*uy,
        dnda*Ds*ux*uy,
        dnda*Ds*(1 + uy^2),
        0,
        (Ds*h*uy)/T,
        0,
        1,
        0,
        (Ds*h*uy)/T))

    source = SVector{3}((
        dxux*n + dyuy*n + (dtux*nux + dtuy*nuy)/utau + dndT*dtT*utau + dndT*dxT*ux + dndT*dyT*uy + (n*(dtux*ux + dtuy*uy))/utau - ((dtux*ux + dtuy*uy)*(nux*ux + nuy*uy))/utau^3 + (nux*ux + nuy*uy + n*(1 + ux^2 + uy^2))/(tau*utau),
        nux + (Ds*(dhdT*nux*(utau*(dxT*ux + dyT*uy) + dtT*(1 + ux^2 + uy^2)) + (n + dndT*T)*(utau*(dxT + dxT*ux^2 + dyT*ux*uy) + dtT*ux*(1 + ux^2 + uy^2))))/(T*utau) + (Ds*h*(nux*(dxuy*tau*utau*ux^3*uy + ux^4*(1 + dtuy*tau*uy) - dyux*tau*utau*ux*uy*(1 + uy^2) - dtux*tau*ux*(-1 + uy^2)*(1 + ux^2 + uy^2) + (1 + uy^2)*(1 + 2*dxux*tau*utau + dyuy*tau*utau + dtuy*tau*uy + uy^2) + ux^2*(2 + dxux*tau*utau + dyuy*tau*utau + (2 - dxux*tau*utau + dyuy*tau*utau)*uy^2 + dtuy*tau*uy*(2 + uy^2))) + tau*(-(dtuy*nuy*ux*(1 + ux^2)*(1 + ux^2 + uy^2)) + dtux*(1 + ux^2 + uy^2)*(nuy*(1 + ux^2)*uy + n*(1 + ux^2 + uy^2)) + utau*(dxux*nuy*ux^3*uy - nuy*ux*(1 + ux^2)*(dxuy*ux + dyuy*uy) + dyux*nuy*(1 + ux^2)*(1 + uy^2) + dxux*n*ux*(1 + ux^2 + uy^2) + dyux*n*uy*(1 + ux^2 + uy^2)))))/(T*tau*utau^3),
        nuy + (Ds*(dhdT*nuy*(utau*(dxT*ux + dyT*uy) + dtT*(1 + ux^2 + uy^2)) + (n + dndT*T)*(dtT*uy*(1 + ux^2 + uy^2) + utau*(dyT + dxT*ux*uy + dyT*uy^2))))/(T*utau) + (Ds*h*(nuy*(1 + dxux*tau*utau + 2*dyuy*tau*utau + dtuy*tau*uy - dxuy*tau*utau*ux^3*uy + 2*uy^2 + dxux*tau*utau*uy^2 + dyuy*tau*utau*uy^2 + dtuy*tau*uy^3 + uy^4 + ux^4*(1 - dtuy*tau*uy) + dtux*tau*ux*(1 + uy^2)*(1 + ux^2 + uy^2) + tau*utau*ux*uy*(-dxuy + dyux*uy^2) + ux^2*(2 + dxux*tau*utau + 2*dyuy*tau*utau + (2 + dxux*tau*utau - dyuy*tau*utau)*uy^2 - dtuy*tau*uy^3)) + tau*(-(dtux*nux*uy*(1 + uy^2)*(1 + ux^2 + uy^2)) + dtuy*(1 + ux^2 + uy^2)*(nux*ux*(1 + uy^2) + n*(1 + ux^2 + uy^2)) + utau*(dyuy*nux*ux*uy^3 + dxuy*nux*(1 + ux^2)*(1 + uy^2) - nux*uy*(dxux*ux + dyux*uy)*(1 + uy^2) + dxuy*n*ux*(1 + ux^2 + uy^2) + dyuy*n*uy*(1 + ux^2 + uy^2)))))/(T*tau*utau^3)))

    return (At, Ax, Ay, source)
end
# end of generated code
