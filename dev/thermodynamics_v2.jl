module FluidumThermodynamicsPrototype

using ForwardDiff
using LinearAlgebra
using StaticArrays

export AbstractEquationOfState,
       GrandCanonicalPoint,
       CanonicalPoint,
       MicrocanonicalPoint,
       ThermodynamicState,
       pressure,
       thermodynamic,
       entropy_density,
       charge_densities,
       charge_density,
       susceptibility,
       susceptibilities,
       energy_density,
       canonical,
       microcanonical,
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

# =============================================================================
# Thermodynamic representations
# =============================================================================

"""
    GrandCanonicalPoint(T, μ)

Grand-canonical coordinates

    (T, μ₁, ..., μₙ).

These coordinates select a point at which the pressure is evaluated. They are
not themselves the complete thermodynamic state.
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

"""
    CanonicalPoint(T, n)

Canonical local thermodynamic coordinates `(T, n₁, ..., nₙ)`.
"""
struct CanonicalPoint{T,N}
    T::T
    n::SVector{N,T}
end

function CanonicalPoint(T::S, n::SVector{N,U}) where {S,U,N}
    R = promote_type(S,U)
    CanonicalPoint{R,N}(R(T), SVector{N,R}(n))
end

"""
    MicrocanonicalPoint(energy_density, n)

Microcanonical local thermodynamic coordinates `(ε, n₁, ..., nₙ)`.

Using energy density here keeps the ensemble terminology literal: the
microcanonical representation fixes energy and conserved charges. Entropy is a
derived thermodynamic quantity.
"""
struct MicrocanonicalPoint{T,N}
    energy_density::T
    n::SVector{N,T}
end

function MicrocanonicalPoint(ε::S, n::SVector{N,U}) where {S,U,N}
    R = promote_type(S,U)
    MicrocanonicalPoint{R,N}(R(ε), SVector{N,R}(n))
end

@inline coordinates(x::GrandCanonicalPoint) = SVector(x.T, x.μ...)

@inline function grandcanonical_from_coordinates(y::SVector{D,T}) where {D,T}
    if D == 1
        return GrandCanonicalPoint(y[1])
    end
    μ = SVector{D-1,T}(ntuple(i -> y[i+1], D-1))
    GrandCanonicalPoint(y[1], μ)
end

# =============================================================================
# Equation of state
# =============================================================================

abstract type AbstractEquationOfState end

"""
    pressure(eos, x::GrandCanonicalPoint)

Fundamental EOS interface. A concrete EOS defines its pressure as a function of
grand-canonical coordinates. Optimized or tabulated EOS implementations may
also specialize `thermodynamic` directly.
"""
function pressure end

# EOS algebra is defined at the model level. ThermodynamicState itself is just
# evaluated data and deliberately has no symbolic/calculus algebra.
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
@inline pressure(eos::ScaledEOS, x::GrandCanonicalPoint) = eos.factor*pressure(eos.eos,x)
@inline pressure(eos::OppositeEOS, x::GrandCanonicalPoint) = -pressure(eos.eos,x)

# =============================================================================
# Evaluated thermodynamic state
# =============================================================================

"""
    ThermodynamicState

Thermodynamic information obtained by evaluating an EOS at a chosen point.
In grand-canonical coordinates

    x = (T, μ₁, ..., μₙ),

it stores

    p(x),  ∇p(x),  ∇²p(x).

The gradient is the Legendre-dual vector

    ∇p = (s, n₁, ..., nₙ),

not another grand-canonical point. The lower-right Hessian block is the charge
susceptibility matrix `χᵢⱼ = ∂nᵢ/∂μⱼ`.
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

Generic second-order ForwardDiff fallback. A concrete/tabulated EOS can
specialize this method while returning the same `ThermodynamicState` interface.
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

Legendre transform

    ε = -p + T s + μ⋅n.
"""
@inline function energy_density(x::GrandCanonicalPoint, th::ThermodynamicState)
    -th.pressure + x.T*entropy_density(th) + dot(x.μ, charge_densities(th))
end

# Representation changes when the state has already been evaluated.
@inline canonical(x::GrandCanonicalPoint, th::ThermodynamicState) =
    CanonicalPoint(x.T, charge_densities(th))

@inline microcanonical(x::GrandCanonicalPoint, th::ThermodynamicState) =
    MicrocanonicalPoint(energy_density(x,th), charge_densities(th))

# Convenience forms that evaluate the EOS once.
@inline function canonical(eos::AbstractEquationOfState, x::GrandCanonicalPoint)
    th = thermodynamic(eos,x)
    canonical(x,th)
end

@inline function microcanonical(eos::AbstractEquationOfState, x::GrandCanonicalPoint)
    th = thermodynamic(eos,x)
    microcanonical(x,th)
end

# =============================================================================
# EOS inversion
# =============================================================================

"""
    grandcanonical(eos, target::CanonicalPoint, μguess; ...)

Invert `(T,n) -> (T,μ)`. At fixed temperature the Newton Jacobian is exactly
the charge-susceptibility matrix `χ`.
"""
function grandcanonical(
    eos::AbstractEquationOfState,
    target::CanonicalPoint{T,N},
    μguess::SVector{N,U};
    atol = 1e-12,
    rtol = 1e-10,
    maxiter::Integer = 30,
) where {T,N,U}
    μ = μguess

    for _ in 1:maxiter
        x = GrandCanonicalPoint(target.T, μ)
        th = thermodynamic(eos,x)
        residual = charge_densities(th) - target.n

        scale = max(norm(target.n,Inf), one(eltype(target.n)))
        if norm(residual,Inf) <= atol + rtol*scale
            return x
        end

        μ -= susceptibilities(th) \ residual
    end

    error("canonical -> grand-canonical inversion did not converge in $maxiter iterations")
end

"""
    grandcanonical(eos, target::MicrocanonicalPoint, guess; ...)

Invert `(ε,n) -> (T,μ)` by Newton iteration.

Only the pressure Hessian is required. Since

    ε = -p + x⋅∇p,       x = (T,μ),

we have

    ∂ε/∂x = (∇²p) x.

The remaining Jacobian rows are the corresponding rows of `∇²p` because
`nᵢ = ∂p/∂μᵢ`.
"""
function grandcanonical(
    eos::AbstractEquationOfState,
    target::MicrocanonicalPoint,
    guess::GrandCanonicalPoint;
    atol = 1e-12,
    rtol = 1e-10,
    maxiter::Integer = 30,
)
    xvec = coordinates(guess)
    D = length(xvec)
    N = D - 1

    length(target.n) == N || throw(DimensionMismatch(
        "guess and microcanonical point have different numbers of charges"
    ))

    target_values = SVector(target.energy_density, target.n...)

    for _ in 1:maxiter
        x = grandcanonical_from_coordinates(xvec)
        th = thermodynamic(eos,x)

        values = SVector(energy_density(x,th), charge_densities(th)...)
        residual = values - target_values

        scale = max(norm(target_values,Inf), one(eltype(target_values)))
        if norm(residual,Inf) <= atol + rtol*scale
            return x
        end

        # Jacobian of (ε,n) with respect to (T,μ).
        dε = th.hessian*xvec
        J = SMatrix{D,D}(ntuple(k -> begin
            i = (k-1) % D + 1
            j = (k-1) ÷ D + 1
            i == 1 ? dε[j] : th.hessian[i,j]
        end, D*D))

        xvec -= J \ residual
    end

    error("microcanonical -> grand-canonical inversion did not converge in $maxiter iterations")
end

# Keep `invert` as a descriptive alias for the inverse representation maps.
@inline invert(eos::AbstractEquationOfState, target::CanonicalPoint, guess; kwargs...) =
    grandcanonical(eos,target,guess; kwargs...)

@inline invert(eos::AbstractEquationOfState, target::MicrocanonicalPoint, guess; kwargs...) =
    grandcanonical(eos,target,guess; kwargs...)

# =============================================================================
# Transport models
# =============================================================================

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
Coupled charge-transport model. Conductivity and relaxation time are matrices,
so off-diagonal transport among conserved charges is supported directly.
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
    eta = model.eta_over_s*entropy_density(th)
    enthalpy = energy_density(x,th) + th.pressure
    tau = eta/(model.Ctau*enthalpy)
    ShearState(eta,tau)
end

@inline function bulk_state(::ZeroBulkTransport, x, th)
    z = zero(th.pressure)
    BulkState(z,z)
end

@inline function bulk_state(model::ConstantZetaOverS, x, th)
    zeta = model.zeta_over_s*entropy_density(th)
    enthalpy = energy_density(x,th) + th.pressure
    tau = zeta/(model.Ctau*enthalpy)
    BulkState(zeta,tau)
end

@inline function charge_state(::ZeroChargeTransport, x::GrandCanonicalPoint{T,N}, th) where {T,N}
    z = zero(th.pressure)
    Z = SMatrix{N,N}(ntuple(_ -> z, N*N))
    ChargeState(Z,Z)
end

@inline charge_state(model::ConstantChargeTransport, x, th) =
    ChargeState(model.kappa,model.tau)

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

# =============================================================================
# Analytic EOS used to exercise the prototype
# =============================================================================

"""
    PolynomialMultiChargeEOS(a, chi)

Toy multi-charge EOS

    p(T,μ) = a T⁴ + 1/2 T² μᵀ chi μ.

It has nontrivial temperature-charge and charge-charge Hessian blocks and is
therefore useful for exercising multi-charge thermodynamics and inversion.
"""
struct PolynomialMultiChargeEOS{A,K} <: AbstractEquationOfState
    a::A
    chi::K
end

@inline function pressure(eos::PolynomialMultiChargeEOS, x::GrandCanonicalPoint)
    length(x.μ) == size(eos.chi,1) || throw(DimensionMismatch(
        "EOS and point have different numbers of charges"
    ))
    eos.a*x.T^4 + (x.T^2/2)*dot(x.μ, eos.chi*x.μ)
end

end # module
