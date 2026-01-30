"""
Conformal Equation of State (EOS) implementation.
"""
Base.@kwdef struct ConformalEOS{T} <: EquationOfState
    g_eff::T = 40.0
end

const π2 = π^2
@inline a_SB(eos::ConformalEOS) = (π2/90) * eos.g_eff

@inline @fastmath pressure(T, eos::ConformalEOS) = a_SB(eos) * T^4 * fmGeV^3

thermodynamic(T::N, eos::ConformalEOS{S}) where {N,S} = begin
    TT = promote_type(N,S)

    a = TT(a_SB(eos)) * TT(fmGeV^3)
    T2 = T*T
    T3 = T2*T
    T4 = T2*T2

    p    = a * T4
    dp   = 4a * T3
    d2p  = 12a * T2

    Thermodynamic{TT,1,1}(p, (dp,), (d2p,))
end