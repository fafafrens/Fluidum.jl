Base.@kwdef struct TanhGEffectiveEOS{T} <: EquationOfState
    g_low::T  = 3.0    # effective low-T dof (tune)
    g_high::T = 40.0   # high-T dof (match your conformal value)
    Tc::T     = 0.17   # GeV, crossover midpoint (tune)
    Δ::T      = 0.03   # GeV, crossover width (tune)
end

@inline a_SB_from_g(g) = (π2/90) * g

@inline @fastmath function g_eff(T, eos::TanhGEffectiveEOS)
    u = (T - eos.Tc) / eos.Δ
    return eos.g_low + (eos.g_high - eos.g_low) * 0.5 * (1 + tanh(u))
end

@inline @fastmath function g_eff_derivative(T, eos::TanhGEffectiveEOS)
    u  = (T - eos.Tc) / eos.Δ
    du = inv(eos.Δ)
    sech2 = inv(cosh(u))^2
    return (eos.g_high - eos.g_low) * 0.5 * sech2 * du
end

@inline @fastmath function g_eff_derivative2(T, eos::TanhGEffectiveEOS)
    u   = (T - eos.Tc) / eos.Δ
    du  = inv(eos.Δ)
    sech2 = inv(cosh(u))^2
    # d/du(sech^2 u) = -2 sech^2 u tanh u
    return (eos.g_high - eos.g_low) * 0.5 * (-2 * sech2 * tanh(u)) * du^2
end

@inline @fastmath function pressure(T, eos::TanhGEffectiveEOS)
    g = g_eff(T, eos)
    a = a_SB_from_g(g)
    return (a * T^4) * fmGeV^3
end

@inline @fastmath function pressure_derivative(T, ::Val{1}, eos::TanhGEffectiveEOS)
    g  = g_eff(T, eos)
    dg = g_eff_derivative(T, eos)
    a0 = (π2/90)
    return (a0 * (dg*T^4 + g*4T^3)) * fmGeV^3
end

@inline @fastmath function pressure_derivative(T, ::Val{2}, eos::TanhGEffectiveEOS)
    g   = g_eff(T, eos)
    dg  = g_eff_derivative(T, eos)
    d2g = g_eff_derivative2(T, eos)
    a0  = (π2/90)
    return (a0 * (d2g*T^4 + 8dg*T^3 + g*12T^2)) * fmGeV^3
end
