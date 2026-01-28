Base.@kwdef struct RunningConformalEOS{T} <: EquationOfState
    g_eff::T = 40.0
    δ::T     = 1e-4     # dimensionless
    T0::T    = 0.2     # GeV (choose a fixed scale)
end

@inline a_SB(eos::RunningConformalEOS) = (π2/90) * eos.g_eff

@inline @fastmath function pressure(T, eos::RunningConformalEOS)
    a = a_SB(eos)
    f = 1 + eos.δ*(eos.T0/T)^2
    return (a*T^4*f) * fmGeV^3
end

@inline @fastmath function pressure_derivative(T, ::Val{1}, eos::RunningConformalEOS)
    a = a_SB(eos)
    δ, T0 = eos.δ, eos.T0
    # P = a(T^4 + δ T0^2 T^2)
    return a*(4T^3 + 2δ*T0^2*T) * fmGeV^3
end

@inline @fastmath function pressure_derivative(T, ::Val{2}, eos::RunningConformalEOS)
    a = a_SB(eos)
    δ, T0 = eos.δ, eos.T0
    return a*(12T^2 + 2δ*T0^2) * fmGeV^3
end
