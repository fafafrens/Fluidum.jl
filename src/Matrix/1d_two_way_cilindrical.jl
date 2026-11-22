# the convention here are T, ur, \[Pi]phiphi, \[Pi]etaeta, \[Pi]B, α, nur, \[Pi]phiphi_Q, \[Pi]etaeta_Q, \[Pi]B_Q (1+1D cylindrical, tau and r)
# with a ConformalEOS there is no \[Pi]B: T, ur, \[Pi]phiphi, \[Pi]etaeta, α, nur, \[Pi]phiphi_Q, \[Pi]etaeta_Q, \[Pi]B_Q
# Nabla_mu (T_m + T_Q)^{mu nu} = 0 with T_Q = (n h - P0) u u + P0 Delta + h (u nu + nu u), P0 = nT, the medium pi and
# Pi rows, and the charm rows of matrix1d_visc_HQ_BG_consistent_m2!; the second moment does not enter T_Q

"""
    matrix1d_two_way!(; species, delta_pipi = true, edge = ColdEdge(), dmn_eps = 1e-8, α_max = 200.0)
The 1+1D two-way system, medium and charm evolved together: layout `two_way_1d` (`two_way_conformal_1d`
for a `ConformalEOS`). `species` is the heavy-particle list (`charm_quark_list()` for one charm quark),
`delta_pipi` as in `matrix1d_visc!`, `edge = nothing` switches the cold-edge decoupling of D_s off.
"""
function matrix1d_two_way!(; species, delta_pipi::Bool = true, edge = ColdEdge(), dmn_eps = 1e-8, α_max = 200.0)
    (A_i, Source, ϕ, tau, X, params) -> matrix1d_two_way!(A_i, Source, ϕ, tau, X, params;
        species = species, delta_pipi = delta_pipi, edge = edge, dmn_eps = dmn_eps, α_max = α_max)
end

function matrix1d_two_way!(A_i, Source, ϕ, tau, X, params; species, delta_pipi::Bool = true, edge = ColdEdge(),
                           dmn_eps = 1e-8, α_max = 200.0)

    T = ϕ[1]
    uτ = sqrt(1 + ϕ[2]^2)

    if params.eos isa ConformalEOS
        ch = charm_eos(T, clamp(ϕ[5], -α_max, α_max), species)
        Ds = DsT(params.diffusion, T) / T / fmGeV * edge_factor(T, edge)

        a = conformal_aSB(params.eos)
        etaVisc = viscosity(T, 4a * T^3, params.shear)
        tauS = τ_shear(T, 4a * T^3, params.shear)
        dpp = delta_pipi ? (4/3) * tauS : zero(tauS)       # delta_pipi(tauS), as in matrix1d_visc!

        (At, Ax, source) = one_d_two_way_matrix_conformal(ϕ, tau, X[1], a, etaVisc, tauS, dpp,
                                                          ch.n, ch.dndT, ch.dnda, ch.h, ch.dhdT, ch.I31, ch.dI31dT, ch.dI31da,
                                                          ch.aQ, ch.bQ, ch.dbQdT, ch.cQ, ch.dcQdT, Ds)
        # vacuum regulator: + dmn_eps u.d α on the number row's own α column only, so no other row and
        # not κ sees it (the background drivers add dmn_eps to dn/dα, which also enters κ = D_s dn/dα)
        At = Base.setindex(At, At[5, 5] + dmn_eps * uτ, 5, 5)
        Ax = Base.setindex(Ax, Ax[5, 5] + dmn_eps * ϕ[2], 5, 5)
    else
        ch = charm_eos(T, clamp(ϕ[6], -α_max, α_max), species)
        Ds = DsT(params.diffusion, T) / T / fmGeV * edge_factor(T, edge)

        therm = thermodynamic(T, params.eos)
        dpt = therm.pressure_derivative[1]
        dptt = therm.pressure_hessian[1]
        etaVisc = viscosity(T, dpt, params.shear)
        tauS = τ_shear(T, dpt, params.shear)
        zeta = bulk_viscosity(T, dpt, params.bulk)
        tauB = τ_bulk(T, dpt, dptt, params.bulk)
        dpp = delta_pipi ? (4/3) * tauS : zero(tauS)

        (At, Ax, source) = one_d_two_way_matrix(ϕ, tau, X[1], therm.pressure, dpt, dptt, zeta, etaVisc, tauS, tauB, dpp,
                                                ch.n, ch.dndT, ch.dnda, ch.h, ch.dhdT, ch.I31, ch.dI31dT, ch.dI31da,
                                                ch.aQ, ch.bQ, ch.dbQdT, ch.cQ, ch.dcQdT, Ds)
        At = Base.setindex(At, At[6, 6] + dmn_eps * uτ, 6, 6)
        Ax = Base.setindex(Ax, Ax[6, 6] + dmn_eps * ϕ[2], 6, 6)
    end

    Ainv = inv(At)
    A_mul_B!(A_i[1], Ainv, Ax)
    jgemvavx!(Source, Ainv, source)
    return nothing
end


# the corona update (ColdEdge(; corona = true)): below T_edge the quasi-linear α update is blended into a
# flux-form update of the charm number, d(τ vol N^τ)/dτ = −(flux differences)/dx with a local Lax–Friedrichs
# flux, which keeps n positive at an under-resolved charm edge; nothing changes for T ≥ T_edge

# n, ∂_T n, ∂_α n of the species list
@inline function charm_density(T, α, species)
    th = thermodynamic(T, α, species)
    return th.pressure, th.pressure_derivative[1], th.pressure_derivative[2]
end

# the LLF flux; a = the larger drift |N^r|/N^τ of the two cells, capped at 1 so a near-vacuum cell is not stiff
@inline function llf_flux(UL, UR, FL, FR)
    a = min(one(UL), max(UL > 0 ? abs(FL) / UL : one(UL), UR > 0 ? abs(FR) / UR : one(UR)))
    return (FL + FR) / 2 - a * (UR - UL) / 2
end

# α̇ that gives the target d(τ vol N^τ)/dτ (vol = r in 1+1D, 1 in 2+1D), given the cell's state and other rates
@inline function corona_alpha_rate(target, vol, τ, T, α, uτ, u̇τ, dντ, Ṫ, Nt, species, dmn_eps, α_max)
    n, dndT, dnda = charm_density(T, clamp(α, -α_max, α_max), species)
    dna = abs(α) < α_max ? dnda : zero(dnda)
    return ((target - vol * Nt) / (τ * vol) - dndT * Ṫ * uτ - n * u̇τ - dντ) / (uτ * (dna + dmn_eps))
end

# N^τ, N^r of 1+1D cell i
@inline function charm_current_1d(u, i, species, iα, α_max)
    ur = u[2, i]; uτ = sqrt(1 + ur^2)
    n, _, _ = charm_density(u[1, i], clamp(u[iα, i], -α_max, α_max), species)
    return n * uτ + ur * u[iα+1, i] / uτ, n * ur + u[iα+1, i]
end

function corona_update_1d!(du, u, τ, rr, dr, iα, species, edge, dmn_eps, α_max)
    iν = iα + 1
    @inbounds for i in 2:size(u, 2)-1
        w = corona_weight(u[1, i], edge)
        iszero(w) && continue
        NtC, NrC = charm_current_1d(u, i, species, iα, α_max)
        rm = rr[i] - dr / 2; rp = rr[i] + dr / 2
        Gm = if rm <= 0
            zero(NtC)
        else
            NtL, NrL = charm_current_1d(u, i - 1, species, iα, α_max); τ * rm * llf_flux(NtL, NtC, NrL, NrC)
        end
        NtR, NrR = charm_current_1d(u, i + 1, species, iα, α_max)
        Gp = τ * rp * llf_flux(NtC, NtR, NrC, NrR)
        ur = u[2, i]; uτ = sqrt(1 + ur^2); ν = u[iν, i]
        u̇ = du[2, i]; ν̇ = du[iν, i]; u̇τ = ur * u̇ / uτ
        dντ = (u̇ * ν + ur * ν̇) / uτ - ur * ν * u̇τ / uτ^2
        αff = corona_alpha_rate(-(Gp - Gm) / dr, rr[i], τ, u[1, i], u[iα, i], uτ, u̇τ, dντ, du[1, i], NtC, species, dmn_eps, α_max)
        du[iα, i] = (1 - w) * du[iα, i] + w * αff
    end
    return nothing
end

"""
    two_way_problem(d, params, init, tspan; species, delta_pipi = true, edge = ColdEdge(), dmn_eps = 1e-8, α_max = 200.0)
The `ODEProblem` of the two-way system on `d` (1+1D or 2+1D; layout `two_way_1d`/`two_way_2d`, or the
conformal ones for a `ConformalEOS`): `problem` with `matrix1d_two_way!`/`matrix2d_two_way!`, plus the
corona update when `edge.corona` (the default). Keywords as in `matrix1d_two_way!`.
"""
function two_way_problem(d::DiscreteFields{T,2,1}, params, init, tspan; species, delta_pipi::Bool = true,
                         edge = ColdEdge(), dmn_eps = 1e-8, α_max = 200.0) where {T}
    base = problem(d, matrix1d_two_way!(; species, delta_pipi, edge, dmn_eps, α_max), params, init, tspan)
    (edge isa ColdEdge && edge.corona) || return base
    iα = params.eos isa ConformalEOS ? 5 : 6
    X = d.discretization
    rr = [X.grid[i][1] for i in eachindex(X.grid)]
    dr = float(Δx(X)[1])
    f = base.f
    function rhs!(du, u, p, t)
        f(du, u, p, t)
        corona_update_1d!(du, u, t, rr, dr, iα, species, edge, dmn_eps, α_max)
        boundary_condition!(du, d)
        return nothing
    end
    return ODEProblem(rhs!, base.u0, base.tspan)
end

# generated by symbolic computation (Mathematica)
"""
    one_d_two_way_matrix(u, tau, rr, pres, dtp, dtdtp, zeta, visc, tauS, tauB, dpp, n, dndT, dnda, h, dhdT, I31, dI31dT, dI31da, aQ, bQ, dbQdT, cQ, dcQdT, Ds)
The two-way system, Nabla_mu (T_m + T_Q)^{mu nu} = 0 with the medium and charm rows; 1+1D cylindrical.
Fields u = (T, ur, piphiphi, pietaeta, piB, α, nur, piphiphi_Q, pietaeta_Q, piB_Q); returns (A_t, A_r, source).
"""
function one_d_two_way_matrix(u, tau, rr, pres, dtp, dtdtp, zeta, visc, tauS, tauB, dpp, n, dndT, dnda, h, dhdT, I31, dI31dT, dI31da, aQ, bQ, dbQdT, cQ, dcQdT, Ds)
    T, ur, piphiphi, pietaeta, piB, α, nur, piphiphi_Q, pietaeta_Q, piB_Q = u[1], u[2], u[3], u[4], u[5], u[6], u[7], u[8], u[9], u[10]
    utau = sqrt(1 + ur^2)
    At = SMatrix{10,10}((
        dtdtp*T + 2*dhdT*nur*ur + dtp*ur^2 + dtdtp*T*ur^2 + n*(-1 + dhdT + dhdT*ur^2) + dndT*(h - T + h*ur^2),
        ((dtp + dndT*h + dtdtp*T)*ur*(1 + ur^2) + dhdT*(nur + 2*nur*ur^2 + n*ur*(1 + ur^2)))/utau,
        0,
        0,
        0,
        dndT*utau,
        (Ds*(dhdT*nur + (n + dndT*T)*ur)*utau)/T,
        (Ds*(-2*dbQdT*nur*ur + 3*dcQdT*piphiphi_Q*(1 + ur^2)))/(6*T*utau),
        (Ds*(-2*dbQdT*nur*ur + 3*dcQdT*pietaeta_Q*(1 + ur^2)))/(6*T*utau),
        (Ds*(5*dbQdT*nur*ur + 3*dI31dT*(1 + ur^2) + 3*dcQdT*piB_Q*(1 + ur^2)))/(6*T*utau),
        2*((piB - pietaeta - piphiphi + dtp*T)*ur + h*(nur + n*ur)),
        ((piB - pietaeta - piphiphi + dtp*T)*(1 + 3*ur^2 + 2*ur^4) + h*(n + 3*nur*ur + 3*n*ur^2 + 2*nur*ur^3 + 2*n*ur^4))/utau^3,
        (ur*(3*dpp*piphiphi - 2*visc))/(3*rr^2*utau),
        (ur*(3*dpp*pietaeta - 2*visc))/(3*tau^2*utau),
        (ur*zeta)/utau,
        (nur + n*ur*(1 + ur^2))/utau^3,
        (Ds*h*(n + nur*ur + n*ur^2))/(T*utau),
        -1/6*(Ds*(2*bQ*nur + 2*aQ*nur*(1 + ur^2) + (2*I31 + 2*cQ*piB_Q - 2*cQ*pietaeta_Q - 5*cQ*piphiphi_Q)*ur*(1 + ur^2)))/(T*utau^3),
        -1/6*(Ds*(2*bQ*nur + 2*aQ*nur*(1 + ur^2) + (2*I31 + 2*cQ*piB_Q - 5*cQ*pietaeta_Q - 2*cQ*piphiphi_Q)*ur*(1 + ur^2)))/(T*utau^3),
        (Ds*(5*bQ*nur + 2*aQ*nur*(1 + ur^2) + (5*I31 + 5*cQ*piB_Q - 2*cQ*(pietaeta_Q + piphiphi_Q))*ur*(1 + ur^2)))/(6*T*utau^3),
        -ur^2,
        -(ur*utau),
        (tauS*utau)/rr^2,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        -ur^2,
        -(ur*utau),
        0,
        (tauS*utau)/tau^2,
        0,
        0,
        0,
        0,
        0,
        0,
        ur^2,
        ur*utau,
        0,
        0,
        tauB*utau,
        0,
        0,
        0,
        0,
        0,
        dnda*(h - T + h*ur^2),
        dnda*h*ur*utau,
        0,
        0,
        0,
        dnda*utau,
        dnda*Ds*ur*utau,
        0,
        0,
        (dI31da*Ds*utau)/(2*T),
        2*h*ur,
        (h*(1 + 2*ur^2))/utau,
        0,
        0,
        0,
        ur/utau,
        (Ds*h*utau)/T,
        -1/3*(bQ*Ds*ur)/(T*utau),
        -1/3*(bQ*Ds*ur)/(T*utau),
        (5*bQ*Ds*ur)/(6*T*utau),
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        (cQ*Ds*utau)/(2*T),
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        (cQ*Ds*utau)/(2*T),
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        (cQ*Ds*utau)/(2*T)))

    Ax = SMatrix{10,10}((
        ((dtp + dndT*h + dtdtp*T)*ur*(1 + ur^2) + dhdT*(nur + 2*nur*ur^2 + n*ur*(1 + ur^2)))/utau,
        dtp + n + dndT*T + 2*dhdT*nur*ur + dtp*ur^2 + dndT*h*ur^2 + dhdT*n*ur^2 + dtdtp*T*ur^2,
        0,
        0,
        0,
        dndT*ur,
        (Ds*(dhdT*nur*ur + n*(1 + ur^2) + dndT*T*(1 + ur^2)))/T,
        -1/6*(Ds*(2*dbQdT*nur - 3*dcQdT*piphiphi_Q*ur))/T,
        -1/6*(Ds*(2*dbQdT*nur - 3*dcQdT*pietaeta_Q*ur))/T,
        (Ds*(5*dbQdT*nur + 3*(dI31dT + dcQdT*piB_Q)*ur))/(6*T),
        ((piB - pietaeta - piphiphi + dtp*T)*(1 + 3*ur^2 + 2*ur^4) + h*(n + 3*nur*ur + 3*n*ur^2 + 2*nur*ur^3 + 2*n*ur^4))/utau^3,
        2*((piB - pietaeta - piphiphi + dtp*T)*ur + h*(nur + n*ur)),
        (3*dpp*piphiphi - 2*visc)/(3*rr^2),
        (3*dpp*pietaeta - 2*visc)/(3*tau^2),
        zeta,
        n,
        (Ds*h*(n*ur*(1 + ur^2) + nur*(2 + ur^2)))/(T*utau^2),
        -1/6*(Ds*(2*aQ*nur*ur + 2*I31*(1 + ur^2) + cQ*(2*piB_Q - 2*pietaeta_Q - 5*piphiphi_Q)*(1 + ur^2)))/(T*utau^2),
        -1/6*(Ds*(2*aQ*nur*ur + 2*I31*(1 + ur^2) + cQ*(2*piB_Q - 5*pietaeta_Q - 2*piphiphi_Q)*(1 + ur^2)))/(T*utau^2),
        (Ds*(2*aQ*nur*ur + 5*I31*(1 + ur^2) + cQ*(5*piB_Q - 2*(pietaeta_Q + piphiphi_Q))*(1 + ur^2)))/(6*T*utau^2),
        -(ur*utau),
        -1 - ur^2,
        (tauS*ur)/rr^2,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        -(ur*utau),
        -1 - ur^2,
        0,
        (tauS*ur)/tau^2,
        0,
        0,
        0,
        0,
        0,
        0,
        ur*utau,
        1 + ur^2,
        0,
        0,
        tauB*ur,
        0,
        0,
        0,
        0,
        0,
        dnda*h*ur*utau,
        dnda*(T + h*ur^2),
        0,
        0,
        0,
        dnda*ur,
        dnda*Ds*(1 + ur^2),
        0,
        0,
        (dI31da*Ds*ur)/(2*T),
        (h*(1 + 2*ur^2))/utau,
        2*h*ur,
        0,
        0,
        0,
        1,
        (Ds*h*ur)/T,
        -1/3*(bQ*Ds)/T,
        -1/3*(bQ*Ds)/T,
        (5*bQ*Ds)/(6*T),
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        (cQ*Ds*ur)/(2*T),
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        (cQ*Ds*ur)/(2*T),
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        (cQ*Ds*ur)/(2*T)))

    source = SVector{10}((
        (-(pietaeta*tau*ur) - piphiphi*tau*ur + dtp*T*tau*ur - pietaeta*tau*ur^3 - piphiphi*tau*ur^3 + dtp*T*tau*ur^3 + pietaeta*rr*utau + dtp*rr*T*utau - pietaeta*rr*ur^2*utau - piphiphi*rr*ur^2*utau + dtp*rr*T*ur^2*utau + h*n*(1 + ur^2)*(tau*ur + rr*utau) + piB*(1 + ur^2)*(tau*ur + rr*utau) + h*nur*(tau + 2*tau*ur^2 + 2*rr*ur*utau))/(rr*tau*utau),
        -((piB + piphiphi + pres)/rr) - (n*T)/rr + (n*T + 2*h*nur*ur + h*n*ur^2)/rr + (piB + pres + (piB + dtp*T)*ur^2 - (pietaeta + piphiphi)*(1 + ur^2))/rr + (h*(nur + 2*nur*ur^2 + n*ur*(1 + ur^2)))/(tau*utau) + ((piB - pietaeta - piphiphi + dtp*T)*ur*utau)/tau,
        (3*piphiphi*(dpp*tau*ur + rr*(tau + dpp*utau)) + 4*tau*ur*visc - 2*rr*utau*visc)/(3*rr^3*tau),
        (3*pietaeta*(dpp*tau*ur + rr*(tau + dpp*utau)) - 2*tau*ur*visc + 4*rr*utau*visc)/(3*rr*tau^3),
        piB + (ur*zeta)/rr + (utau*zeta)/tau,
        (nur + n*ur)/rr + (n + nur*ur + n*ur^2)/(tau*utau),
        nur + (Ds*h*nur*(rr + rr*ur^2 + tau*ur*utau))/(rr*T*tau*utau),
        piphiphi_Q + (Ds*(-2*bQ*nur*(rr*ur - 2*tau*utau) - 2*I31*(rr + rr*ur^2 - 2*tau*ur*utau) + cQ*(-2*pietaeta_Q*rr*(1 + ur^2) + 3*piphiphi_Q*rr*(1 + ur^2) + 7*piphiphi_Q*tau*ur*utau - 2*piB_Q*(rr + rr*ur^2 - 2*tau*ur*utau))))/(6*rr*T*tau*utau),
        pietaeta_Q + (Ds*(4*I31*rr*(1 + ur^2) - 2*I31*tau*ur*utau + 2*bQ*nur*(2*rr*ur - tau*utau) + cQ*(4*piB_Q*rr*(1 + ur^2) + 7*pietaeta_Q*rr*(1 + ur^2) - 2*piB_Q*tau*ur*utau + 3*pietaeta_Q*tau*ur*utau - 2*piphiphi_Q*tau*ur*utau)))/(6*rr*T*tau*utau),
        piB_Q + (Ds*(5*bQ*nur*(rr*ur + tau*utau) + 5*I31*(rr + rr*ur^2 + tau*ur*utau) + cQ*(2*pietaeta_Q*rr*(1 + ur^2) + 2*piphiphi_Q*tau*ur*utau + 5*piB_Q*(rr + rr*ur^2 + tau*ur*utau))))/(6*rr*T*tau*utau)))

    return (At, Ax, source)
end


"""
    one_d_two_way_matrix_conformal(u, tau, rr, aSB, visc, tauS, dpp, n, dndT, dnda, h, dhdT, I31, dI31dT, dI31da, aQ, bQ, dbQdT, cQ, dcQdT, Ds)
The two-way system, Nabla_mu (T_m + T_Q)^{mu nu} = 0 with the medium and charm rows; 1+1D cylindrical, conformal medium.
Fields u = (T, ur, piphiphi, pietaeta, α, nur, piphiphi_Q, pietaeta_Q, piB_Q); returns (A_t, A_r, source).
"""
function one_d_two_way_matrix_conformal(u, tau, rr, aSB, visc, tauS, dpp, n, dndT, dnda, h, dhdT, I31, dI31dT, dI31da, aQ, bQ, dbQdT, cQ, dcQdT, Ds)
    T, ur, piphiphi, pietaeta, α, nur, piphiphi_Q, pietaeta_Q, piB_Q = u[1], u[2], u[3], u[4], u[5], u[6], u[7], u[8], u[9]
    utau = sqrt(1 + ur^2)
    At = SMatrix{9,9}((
        2*dhdT*nur*ur + 4*aSB*T^3*(3 + 4*ur^2) + n*(-1 + dhdT + dhdT*ur^2) + dndT*(h - T + h*ur^2),
        ((dndT*h + 16*aSB*T^3)*ur*(1 + ur^2) + dhdT*(nur + 2*nur*ur^2 + n*ur*(1 + ur^2)))/utau,
        0,
        0,
        dndT*utau,
        (Ds*(dhdT*nur + (n + dndT*T)*ur)*utau)/T,
        (Ds*(-2*dbQdT*nur*ur + 3*dcQdT*piphiphi_Q*(1 + ur^2)))/(6*T*utau),
        (Ds*(-2*dbQdT*nur*ur + 3*dcQdT*pietaeta_Q*(1 + ur^2)))/(6*T*utau),
        (Ds*(5*dbQdT*nur*ur + 3*dI31dT*(1 + ur^2) + 3*dcQdT*piB_Q*(1 + ur^2)))/(6*T*utau),
        -2*(pietaeta + piphiphi - 4*aSB*T^4)*ur + 2*h*(nur + n*ur),
        (-((pietaeta + piphiphi - 4*aSB*T^4)*(1 + 3*ur^2 + 2*ur^4)) + h*(n + 3*nur*ur + 3*n*ur^2 + 2*nur*ur^3 + 2*n*ur^4))/utau^3,
        (ur*(3*dpp*piphiphi - 2*visc))/(3*rr^2*utau),
        (ur*(3*dpp*pietaeta - 2*visc))/(3*tau^2*utau),
        (nur + n*ur*(1 + ur^2))/utau^3,
        (Ds*h*(n + nur*ur + n*ur^2))/(T*utau),
        -1/6*(Ds*(2*bQ*nur + 2*aQ*nur*(1 + ur^2) + (2*I31 + 2*cQ*piB_Q - 2*cQ*pietaeta_Q - 5*cQ*piphiphi_Q)*ur*(1 + ur^2)))/(T*utau^3),
        -1/6*(Ds*(2*bQ*nur + 2*aQ*nur*(1 + ur^2) + (2*I31 + 2*cQ*piB_Q - 5*cQ*pietaeta_Q - 2*cQ*piphiphi_Q)*ur*(1 + ur^2)))/(T*utau^3),
        (Ds*(5*bQ*nur + 2*aQ*nur*(1 + ur^2) + (5*I31 + 5*cQ*piB_Q - 2*cQ*(pietaeta_Q + piphiphi_Q))*ur*(1 + ur^2)))/(6*T*utau^3),
        -ur^2,
        -(ur*utau),
        (tauS*utau)/rr^2,
        0,
        0,
        0,
        0,
        0,
        0,
        -ur^2,
        -(ur*utau),
        0,
        (tauS*utau)/tau^2,
        0,
        0,
        0,
        0,
        0,
        dnda*(h - T + h*ur^2),
        dnda*h*ur*utau,
        0,
        0,
        dnda*utau,
        dnda*Ds*ur*utau,
        0,
        0,
        (dI31da*Ds*utau)/(2*T),
        2*h*ur,
        (h*(1 + 2*ur^2))/utau,
        0,
        0,
        ur/utau,
        (Ds*h*utau)/T,
        -1/3*(bQ*Ds*ur)/(T*utau),
        -1/3*(bQ*Ds*ur)/(T*utau),
        (5*bQ*Ds*ur)/(6*T*utau),
        0,
        0,
        0,
        0,
        0,
        0,
        (cQ*Ds*utau)/(2*T),
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        (cQ*Ds*utau)/(2*T),
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        (cQ*Ds*utau)/(2*T)))

    Ax = SMatrix{9,9}((
        ((dndT*h + 16*aSB*T^3)*ur*(1 + ur^2) + dhdT*(nur + 2*nur*ur^2 + n*ur*(1 + ur^2)))/utau,
        n + dndT*T + 4*aSB*T^3 + 2*dhdT*nur*ur + dndT*h*ur^2 + dhdT*n*ur^2 + 16*aSB*T^3*ur^2,
        0,
        0,
        dndT*ur,
        (Ds*(dhdT*nur*ur + n*(1 + ur^2) + dndT*T*(1 + ur^2)))/T,
        -1/6*(Ds*(2*dbQdT*nur - 3*dcQdT*piphiphi_Q*ur))/T,
        -1/6*(Ds*(2*dbQdT*nur - 3*dcQdT*pietaeta_Q*ur))/T,
        (Ds*(5*dbQdT*nur + 3*(dI31dT + dcQdT*piB_Q)*ur))/(6*T),
        (-((pietaeta + piphiphi - 4*aSB*T^4)*(1 + 3*ur^2 + 2*ur^4)) + h*(n + 3*nur*ur + 3*n*ur^2 + 2*nur*ur^3 + 2*n*ur^4))/utau^3,
        -2*(pietaeta + piphiphi - 4*aSB*T^4)*ur + 2*h*(nur + n*ur),
        (3*dpp*piphiphi - 2*visc)/(3*rr^2),
        (3*dpp*pietaeta - 2*visc)/(3*tau^2),
        n,
        (Ds*h*(n*ur*(1 + ur^2) + nur*(2 + ur^2)))/(T*utau^2),
        -1/6*(Ds*(2*aQ*nur*ur + 2*I31*(1 + ur^2) + cQ*(2*piB_Q - 2*pietaeta_Q - 5*piphiphi_Q)*(1 + ur^2)))/(T*utau^2),
        -1/6*(Ds*(2*aQ*nur*ur + 2*I31*(1 + ur^2) + cQ*(2*piB_Q - 5*pietaeta_Q - 2*piphiphi_Q)*(1 + ur^2)))/(T*utau^2),
        (Ds*(2*aQ*nur*ur + 5*I31*(1 + ur^2) + cQ*(5*piB_Q - 2*(pietaeta_Q + piphiphi_Q))*(1 + ur^2)))/(6*T*utau^2),
        -(ur*utau),
        -1 - ur^2,
        (tauS*ur)/rr^2,
        0,
        0,
        0,
        0,
        0,
        0,
        -(ur*utau),
        -1 - ur^2,
        0,
        (tauS*ur)/tau^2,
        0,
        0,
        0,
        0,
        0,
        dnda*h*ur*utau,
        dnda*(T + h*ur^2),
        0,
        0,
        dnda*ur,
        dnda*Ds*(1 + ur^2),
        0,
        0,
        (dI31da*Ds*ur)/(2*T),
        (h*(1 + 2*ur^2))/utau,
        2*h*ur,
        0,
        0,
        1,
        (Ds*h*ur)/T,
        -1/3*(bQ*Ds)/T,
        -1/3*(bQ*Ds)/T,
        (5*bQ*Ds)/(6*T),
        0,
        0,
        0,
        0,
        0,
        0,
        (cQ*Ds*ur)/(2*T),
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        (cQ*Ds*ur)/(2*T),
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        (cQ*Ds*ur)/(2*T)))

    source = SVector{9}((
        (-(piphiphi*tau*ur) + 4*aSB*T^4*tau*ur - piphiphi*tau*ur^3 + 4*aSB*T^4*tau*ur^3 + 4*aSB*rr*T^4*utau - piphiphi*rr*ur^2*utau + 4*aSB*rr*T^4*ur^2*utau - pietaeta*(tau*(ur + ur^3) + rr*(-1 + ur^2)*utau) + h*(n*(1 + ur^2)*(tau*ur + rr*utau) + nur*(tau + 2*tau*ur^2 + 2*rr*ur*utau)))/(rr*tau*utau),
        (-(piphiphi*rr*ur) + 4*aSB*rr*T^4*ur - piphiphi*rr*ur^3 + 4*aSB*rr*T^4*ur^3 - 2*piphiphi*tau*utau - piphiphi*tau*ur^2*utau + 4*aSB*T^4*tau*ur^2*utau - pietaeta*(1 + ur^2)*(rr*ur + tau*utau) + h*n*ur*(rr + rr*ur^2 + tau*ur*utau) + h*nur*(rr + 2*rr*ur^2 + 2*tau*ur*utau))/(rr*tau*utau),
        (3*piphiphi*(dpp*tau*ur + rr*(tau + dpp*utau)) + 4*tau*ur*visc - 2*rr*utau*visc)/(3*rr^3*tau),
        (3*pietaeta*(dpp*tau*ur + rr*(tau + dpp*utau)) - 2*tau*ur*visc + 4*rr*utau*visc)/(3*rr*tau^3),
        (nur + n*ur)/rr + (n + nur*ur + n*ur^2)/(tau*utau),
        nur + (Ds*h*nur*(rr + rr*ur^2 + tau*ur*utau))/(rr*T*tau*utau),
        piphiphi_Q + (Ds*(-2*bQ*nur*(rr*ur - 2*tau*utau) - 2*I31*(rr + rr*ur^2 - 2*tau*ur*utau) + cQ*(-2*pietaeta_Q*rr*(1 + ur^2) + 3*piphiphi_Q*rr*(1 + ur^2) + 7*piphiphi_Q*tau*ur*utau - 2*piB_Q*(rr + rr*ur^2 - 2*tau*ur*utau))))/(6*rr*T*tau*utau),
        pietaeta_Q + (Ds*(4*I31*rr*(1 + ur^2) - 2*I31*tau*ur*utau + 2*bQ*nur*(2*rr*ur - tau*utau) + cQ*(4*piB_Q*rr*(1 + ur^2) + 7*pietaeta_Q*rr*(1 + ur^2) - 2*piB_Q*tau*ur*utau + 3*pietaeta_Q*tau*ur*utau - 2*piphiphi_Q*tau*ur*utau)))/(6*rr*T*tau*utau),
        piB_Q + (Ds*(5*bQ*nur*(rr*ur + tau*utau) + 5*I31*(rr + rr*ur^2 + tau*ur*utau) + cQ*(2*pietaeta_Q*rr*(1 + ur^2) + 2*piphiphi_Q*tau*ur*utau + 5*piB_Q*(rr + rr*ur^2 + tau*ur*utau))))/(6*rr*T*tau*utau)))

    return (At, Ax, source)
end
# end of generated code
