module FluidumThermodynamicsPrototype

using ForwardDiff
using LinearAlgebra
using StaticArrays

export AbstractEquationOfState,
       GrandCanonicalPoint,
       ThermodynamicState,
       CanonicalPoint,
       ConservedPoint,
       pressure,
       thermodynamic,
       entropy_density,
       charge_densities,
       charge_density,
       susceptibility,
       susceptibilities,
       energy_density,
       grandcanonical,
       invert,
       AbstractShearTransport,
       AbstractBulkTransport,
       AbstractChargeTransport,
       ConstantEtaOverS,
       ConstantZetaOverS,
       ConstantChargeTransport,
       ZeroShearTransport,
       ZeroBulkTransport,
       ZeroChargeTransport,
       TransportModel,
       FluidModel,
       ShearState,
       BulkState,
       ChargeState,
       TransportState,
       transport,
       PolynomialMultiChargeEOS

# -----------------------------------------------------------------------------
# Thermodynamic coordinates
# -----------------------------------------------------------------------------

"""
    GrandCanonicalPoint(T, μ)

A point in grand-canonical intensive-variable space

    x = (T, μ₁, ..., μₙ).

The chemical potentials are stored in an `SVector`, so the number of conserved
charges is part of the concrete type.
"""
struct GrandCanonicalPoint{T,N}
    T::T
    μ::SVector{N,T}
end

GrandCanonicalPoint(T::R) where {R} =
    GrandCanonicalPoint{R,0}(T, SVector{0,R}())

function GrandCanonicalPoint(T::S, μ::SVector{N,U}) where {S,U,N}
    R = promote_type(S,U)
    GrandCanonicalPoint{R,N}(R(T), SVector{N,R}(μ))
end

GrandCanonicalPoint(T, μ...) = GrandCanonicalPoint(T, SVector(μ...))

@inline coordinates(x::GrandCanonicalPoint) = SVector(x.T, x.μ...)

@inline function grandcanonical_from_coordinates(y::SVector{D,T}) where {D,T}
    if D == 1
        return GrandCanonicalPoint(y[1])
    end
    μ = SVector{D-1,T}(ntuple(i -> y[i+1], D-1))
    GrandCanonicalPoint(y[1], μ)
end

# -----------------------------------------------------------------------------
# EOS interface
# -----------------------------------------------------------------------------

abstract type AbstractEquationOfState end

"""
    pressure(eos, x::GrandCanonicalPoint)

Fundamental EOS interface. Concrete EOS models only need to define the pressure
as a function of the grand-canonical coordinates. An optimized EOS may also
override `thermodynamic` directly.
"""
function pressure end

# EOS-level algebra. There is deliberately no algebra on ThermodynamicState.
struct SumEOS{A<:AbstractEquationOfState,B<:AbstractEquationOfState} <: AbstractEquationOfState
    a::A
    b::B
end

struct ScaledEOS{T,E<:AbstractEquationOfState} <: AbstractEquationOfState
    factor::T
    eos::E
end

struct OppositeEOS{E<:AbstractEquationOfState} <: AbstractEquationOfState
    eos::E
end

Base.:+(a::AbstractEquationOfState, b::AbstractEquationOfState) = SumEOS(a,b)
Base.:-(a::AbstractEquationOfState) = OppositeEOS(a)
Base.:-(a::AbstractEquationOfState, b::AbstractEquationOfState) = a + (-b)
Base.:*(a::Number, eos::AbstractEquationOfState) = ScaledEOS(a,eos)
Base.:*(eos::AbstractEquationOfState, a::Number) = a*eos
Base.:/(eos::AbstractEquationOfState, a::Number) = inv(a)*eos

@inline pressure(eos::SumEOS, x::GrandCanonicalPoint) = pressure(eos.a,x) + pressure(eos.b,x)
@inline pressure(eos::ScaledEOS, x::GrandCanonicalPoint) = eos.factor * pressure(eos.eos,x)
@inline pressure(eos::OppositeEOS, x::GrandCanonicalPoint) = -pressure(eos.eos,x)

# -----------------------------------------------------------------------------
# Evaluated thermodynamics
# -----------------------------------------------------------------------------

"""
    ThermodynamicState

Value, gradient and Hessian of the pressure at a grand-canonical point.
For x = (T, μ₁, ..., μₙ),

    gradient = (s, n₁, ..., nₙ)

and the lower-right Hessian block is the charge-susceptibility matrix.
"""
struct ThermodynamicState{T,D,G,H}
    pressure::T
    gradient::G
    hessian::H
end

function ThermodynamicState(p, g::SVector{D,T}, H::SMatrix{D,D,T}) where {D,T}
    ThermodynamicState{T,D,typeof(g),typeof(H)}(p,g,H)
end

"""
    thermodynamic(eos, x)

Generic second-order ForwardDiff fallback. Concrete/tabulated EOS implementations
can specialize this method and return the same `ThermodynamicState` interface.
"""
function thermodynamic(eos::AbstractEquationOfState, x::GrandCanonicalPoint)
    y = coordinates(x)
    D = length(y)
    f(z) = pressure(eos, grandcanonical_from_coordinates(z))

    p = f(y)
    g = SVector{D}(ForwardDiff.gradient(f,y))
    H = SMatrix{D,D}(ForwardDiff.hessian(f,y))

    ThermodynamicState(p,g,H)
end

@inline entropy_density(th::ThermodynamicState) = th.gradient[1]

@inline function charge_densities(th::ThermodynamicState{T,D}) where {T,D}
    N = D - 1
    SVector{N}(ntuple(i -> th.gradient[i+1], N))
end

@inline charge_density(th::ThermodynamicState, i::Integer) = th.gradient[i+1]
@inline susceptibility(th::ThermodynamicState, i::Integer, j::Integer) = th.hessian[i+1,j+1]

@inline function susceptibilities(th::ThermodynamicState{T,D}) where {T,D}
    N = D - 1
    SMatrix{N,N}(ntuple(k -> begin
        i = (k-1) % N + 1
        j = (k-1) ÷ N + 1
        th.hessian[i+1,j+1]
    end, N*N))
end

"""
    energy_density(x, th)

Full Legendre transform of the pressure,

    ε = -p + T s + μ⋅n.
"""
@inline function energy_density(x::GrandCanonicalPoint, th::ThermodynamicState)
    -th.pressure + x.T*entropy_density(th) + dot(x.μ, charge_densities(th))
end

# -----------------------------------------------------------------------------
# Alternative thermodynamic coordinates and inversion
# -----------------------------------------------------------------------------

"""Coordinates (T, n₁, ..., nₙ), useful for fixed-temperature inversion."""
struct CanonicalPoint{T,N}
    T::T
    n::SVector{N,T}
end

function CanonicalPoint(T::S, n::SVector{N,U}) where {S,U,N}
    R = promote_type(S,U)
    CanonicalPoint{R,N}(R(T), SVector{N,R}(n))
end

"""Coordinates (s, n₁, ..., nₙ), conjugate to (T, μ₁, ..., μₙ)."""
struct ConservedPoint{T,N}
    s::T
    n::SVector{N,T}
end

function ConservedPoint(s::S, n::SVector{N,U}) where {S,U,N}
    R = promote_type(S,U)
    ConservedPoint{R,N}(R(s), SVector{N,R}(n))
end

@inline coordinates(y::ConservedPoint) = SVector(y.s, y.n...)

"""
    grandcanonical(eos, target::CanonicalPoint, μguess; ...)

Solve n(T,μ) = target.n at fixed temperature. The Newton Jacobian is exactly the
susceptibility matrix χᵢⱼ.
"""
function grandcanonical(
    eos::AbstractEquationOfState,
    target::CanonicalPoint{T,N},
    μguess::SVector{N};
    atol = 1e-12,
    rtol = 1e-10,
    maxiter::Integer = 30,
) where {T,N}
    μ = μguess

    for _ in 1:maxiter
        x = GrandCanonicalPoint(target.T, μ)
        th = thermodynamic(eos,x)
        n = charge_densities(th)
        residual = n - target.n

        if norm(residual, Inf) <= atol + rtol*max(norm(target.n,Inf), one(eltype(target.n)))
            return x
        end

        χ = susceptibilities(th)
        μ -= χ \ residual
    end

    error("fixed-T EOS inversion did not converge in $maxiter iterations")
end

"""
    invert(eos, target::ConservedPoint, guess::GrandCanonicalPoint; ...)

Solve ∇p(T,μ) = (s,n) using Newton iteration. The pressure Hessian is the exact
Jacobian of this map.
"""
function invert(
    eos::AbstractEquationOfState,
    target::ConservedPoint,
    guess::GrandCanonicalPoint;
    atol = 1e-12,
    rtol = 1e-10,
    maxiter::Integer = 30,
)
    xvec = coordinates(guess)
    ytarget = coordinates(target)

    length(xvec) == length(ytarget) || throw(DimensionMismatch("guess and target have different thermodynamic dimensions"))

    for _ in 1:maxiter
        x = grandcanonical_from_coordinates(xvec)
        th = thermodynamic(eos,x)
        residual = th.gradient - ytarget

        if norm(residual, Inf) <= atol + rtol*max(norm(ytarget,Inf), one(eltype(ytarget)))
            return x
        end

        xvec -= th.hessian \ residual
    end

    error("full EOS inversion did not converge in $maxiter iterations")
end

# -----------------------------------------------------------------------------
# Transport models
# -----------------------------------------------------------------------------

abstract type AbstractShearTransport end
abstract type AbstractBulkTransport end
abstract type AbstractChargeTransport end

struct ZeroShearTransport <: AbstractShearTransport end
struct ZeroBulkTransport <: AbstractBulkTransport end
struct ZeroChargeTransport <: AbstractChargeTransport end

"""Simple η/s model with τπ = η/[Cτ(ε+p)]."""
struct ConstantEtaOverS{T} <: AbstractShearTransport
    eta_over_s::T
    Ctau::T
end

"""Simple ζ/s model with τΠ = ζ/[Cτ(ε+p)]."""
struct ConstantZetaOverS{T} <: AbstractBulkTransport
    zeta_over_s::T
    Ctau::T
end

"""
A generic coupled charge sector. Both the conductivity and relaxation time are
matrices, so off-diagonal charge transport is supported from the start.
"""
struct ConstantChargeTransport{K,R} <: AbstractChargeTransport
    kappa::K
    tau::R
end

struct TransportModel{S<:AbstractShearTransport,B<:AbstractBulkTransport,C<:AbstractChargeTransport}
    shear::S
    bulk::B
    charge::C
end

struct FluidModel{E<:AbstractEquationOfState,T<:TransportModel}
    eos::E
    transport::T
end

struct ShearState{T}
    eta::T
    tau::T
end

struct BulkState{T}
    zeta::T
    tau::T
end

struct ChargeState{K,R}
    kappa::K
    tau::R
end

struct TransportState{S,B,C}
    shear::S
    bulk::B
    charge::C
end

@inline function shear_state(::ZeroShearTransport, x, th)
    z = zero(th.pressure)
    ShearState(z,z)
end

@inline function shear_state(model::ConstantEtaOverS, x, th)
    eta = model.eta_over_s * entropy_density(th)
    enthalpy = energy_density(x,th) + th.pressure
    tau = eta/(model.Ctau*enthalpy)
    ShearState(eta,tau)
end

@inline function bulk_state(::ZeroBulkTransport, x, th)
    z = zero(th.pressure)
    BulkState(z,z)
end

@inline function bulk_state(model::ConstantZetaOverS, x, th)
    zeta = model.zeta_over_s * entropy_density(th)
    enthalpy = energy_density(x,th) + th.pressure
    tau = zeta/(model.Ctau*enthalpy)
    BulkState(zeta,tau)
end

@inline function charge_state(::ZeroChargeTransport, x::GrandCanonicalPoint{T,N}, th) where {T,N}
    z = zero(th.pressure)
    Z = SMatrix{N,N}(ntuple(_ -> z, N*N))
    ChargeState(Z,Z)
end

@inline charge_state(model::ConstantChargeTransport, x, th) = ChargeState(model.kappa, model.tau)

@inline function transport(model::TransportModel, x::GrandCanonicalPoint, th::ThermodynamicState)
    TransportState(
        shear_state(model.shear,x,th),
        bulk_state(model.bulk,x,th),
        charge_state(model.charge,x,th),
    )
end

@inline function transport(fluid::FluidModel, x::GrandCanonicalPoint)
    th = thermodynamic(fluid.eos,x)
    th, transport(fluid.transport,x,th)
end

# -----------------------------------------------------------------------------
# Small analytic EOS used to exercise the prototype
# -----------------------------------------------------------------------------

"""
    PolynomialMultiChargeEOS(a, chi)

Toy multi-charge EOS

    p(T,μ) = a T⁴ + 1/2 T² μᵀ chi μ.

It is intentionally simple but has nontrivial T-μ and charge-charge Hessian
blocks, making it useful for prototyping multi-charge thermodynamics.
"""
struct PolynomialMultiChargeEOS{A,K} <: AbstractEquationOfState
    a::A
    chi::K
end

@inline function pressure(eos::PolynomialMultiChargeEOS, x::GrandCanonicalPoint)
    length(x.μ) == size(eos.chi,1) || throw(DimensionMismatch("EOS and point have different numbers of charges"))
    eos.a*x.T^4 + (x.T^2/2)*dot(x.μ, eos.chi*x.μ)
end

end # module
