Base.@kwdef struct ConformalEOS{T} <: EquationOfState
    g_eff::T = 40.0    # effective degeneracy
end

const π2 = π^2

@inline a_SB(eos::ConformalEOS) = (π2/90) * eos.g_eff

@inline @fastmath function pressure(T, eos::ConformalEOS)
    a = a_SB(eos)
    return a * T^4 * fmGeV^3
end

@inline @fastmath function pressure_derivative(T, ::Val{1}, eos::ConformalEOS)
    a = a_SB(eos)
    return 4a * T^3 * fmGeV^3
end

@inline @fastmath function pressure_derivative(T, ::Val{2}, eos::ConformalEOS)
    a = a_SB(eos)
    return 12a * T^2 * fmGeV^3
end
