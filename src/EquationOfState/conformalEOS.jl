"""
    ConformalEOS(; g_eff = 40.0)
Conformal medium, P = a_SB T^4 with a_SB = g_eff π²/90 (T in GeV, P in GeV/fm³).
"""
Base.@kwdef struct ConformalEOS{T} <: EquationOfState
    g_eff::T = 40.0
end

const π2 = π^2
@inline a_SB(eos::ConformalEOS) = (π2/90) * eos.g_eff

# P = conformal_aSB(eos) T^4 in GeV/fm³; the aSB argument of the conformal two-way functions
@inline conformal_aSB(eos::ConformalEOS) = a_SB(eos) * fmGeV^3

@inline @fastmath pressure(T, eos::ConformalEOS) = a_SB(eos) * T^4 * fmGeV^3

thermodynamic(T::N, eos::ConformalEOS{S}) where {N,S} = begin
    # float(): an integer T would otherwise give TT = Int
    TT = float(promote_type(N,S))

    a = TT(a_SB(eos)) * TT(fmGeV^3)
    T2 = T*T
    T3 = T2*T
    T4 = T2*T2

    p    = a * T4
    dp   = 4a * T3
    d2p  = 12a * T2

    Thermodynamic{TT,1,1}(p, (dp,), (d2p,))
end