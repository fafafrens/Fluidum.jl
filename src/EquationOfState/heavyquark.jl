Base.@kwdef struct Heavy_Quark{T,S<:CharmHadronList} <:EquationOfState
    mass::T = 1.5
    eos::FluiduMEoS{T} = FluiduMEoS()
    hadron_list::S = HadronResonaceGas_ccbar(readresonancelist(), 30.)
end

function Heavy_Quark(particle_list, ccbar)
    list = HadronResonaceGas_ccbar(readresonancelist(;name_file = particle_list), ccbar)
    return Heavy_Quark(hadron_list = list)
end

@inline @fastmath thermodynamic(T, x::Heavy_Quark) = thermodynamic(T,x.eos)

function free_charm(T,fug,x::Heavy_Quark)
    m = x.mass
    b2 = besselkx(2,m/T)
    b1 = besselk1x(m/T)
    b3 = b1+4/(m/T)*b2

    ex = exp(fug-m/T)
    deg = 6 #spin x2,color x3
    density = deg*(m^2* T /(2 *π^2)* ex* b2)* fmGeV^3; #fm-3
    densityDerT = deg*((ex*m^2*(m*b1 +2*T*b2 + m*b3))/(4*π^2*T))* fmGeV^3 ; #(*fm^-3/GeV*)
    densityDerμ= deg*(m^2* T /(2 *π^2)* ex* b2)* fmGeV^3+0.0001; #(*fm^-3/GeV*)
    densityDerTDerμ =deg*(ex*m^2*(m*b1 + 2*T*b2 + m*b3)/(4*π^2*T))* fmGeV^3 ;
    return (density,densityDerT,densityDerμ,densityDerTDerμ )
end

function free_hadron(T,μ,deg,q; m = 1.5)
    b1 = besselk1x(m/T)
    b2 = besselkx(2,m/T) 
    b3 = b1+4/(m/T)*b2

    ex = exp(q*μ - m/T)  
    density = deg*(T /(2 *π^2)*m^2* ex* b2)* fmGeV3; #fm-3
    densityDerT = deg*((ex*m^2*(m*b1 +2*T*b2 + m*b3))/(4*π^2*T))* fmGeV3 ; #(*fm^-3/GeV)
    
    densityDerμ= density/T+0.0001;
    densityDerTDerμ =deg*(ex*m^2*(m*b1 + 2*T*b2 + m*b3)/(4*π^2*T))* fmGeV3 ;
    return (density,densityDerT,densityDerμ,densityDerTDerμ)
end

@inline function DsT(m::ConstDiffusion,T) 
    m.DsT
end

@inline function DsT(m::LinearDiffusion,T; Tc = m.Tc)
    ifelse(T > Tc, m.slope*T + m.offset, m.slope*Tc + m.offset)
end 


function diffusion(T,n,y::Diffusion)
    density = n #fm-^3
    κ = DsT(y,T)/T*density/fmGeV #fm^-2
end

function diffusion_hadron(T,μ,x::Heavy_Quark,y::Diffusion)
    κ = DsT(y,T)/T*normalization(T,μ,x)/fmGeV
end


function τ_diffusion(T,y::Diffusion) 
    m = y.mass
    b2 = besselkx(2,m/T)*exp(-m/T)
    b1 = besselkx(1,m/T)*exp(-m/T)
   
    b3 = b1+4/(m/T)*b2  
    b4 = b2 + 6/(m/T)*b3  
    b5 = b3+8/(m/T)*b4

  return ((2*π*DsT(y,T)) *m^3/T^3/(96*π*T)*(2*b1 - 3*b3 +b5)/b2)/fmGeV; #
  
end


function τ_diffusion_hadron(T,μ,x::Heavy_Quark,y::Diffusion) 
    tauq = 0
    for i in x.hadron_list.particle_list
        m = i.Mass
        q = i.Nc + i.Nac
        deg = i.Degeneracy
        b2 = besselkx(2,m/T)*exp(-m/T)
        b1 = besselk1x(m/T)*exp(-m/T)
        b3 = b1+4/(m/T)*b2  
        b4 = b2 + 6/(m/T)*b3  
        b5 = b3+8/(m/T)*b4
        ex = exp(q*μ)
        tauq += deg*((2*π*DsT(y,T))/(192*π^3*T^3)*q^2*m^5*ex*(2*b1 - 3*b3 +b5)); 
    end   
  return tauq/normalization(T,μ,x)*(fmGeV^2); #
end


struct QGPViscosity{T}<:ShearViscosity
    ηs::T
    Cs::T
end


@inline function bulk_viscosity(T,entropy,y::SimpleBulkViscosity{N}) where {N}
    y.ζs/(1+((T-0.175)/0.024)^2)*invfmGeV*entropy
 end

@inline function τ_bulk(T,entropy,dtdtp,y::SimpleBulkViscosity{N}) where {N}
     cs2= entropy/(T*dtdtp)
     bulk_viscosity(T,entropy,y)/(T*entropy*y.Cζ)*1/(1/3-cs2)^2+0.1
 end


 @inline function bulk_viscosity(T,y::ZeroBulkViscosity)  
    zero(promote_type(typeof(T)))
end

@inline function bulk_viscosity(T,entropy,y::ZeroBulkViscosity)  
    zero(promote_type(typeof(T)))
end

@inline function τ_bulk(T,y::ZeroBulkViscosity)  
    one(promote_type(typeof(T)))
end

@inline function τ_bulk(T,entropy,y::ZeroBulkViscosity)  
    one(promote_type(typeof(T)))
end

@inline function τ_bulk(T,entropy,dtdtp,y::ZeroBulkViscosity)  
    one(promote_type(typeof(T)))
end 


@inline function viscosity(T,entropy,y::QGPViscosity{N}) where {N}
    y.ηs*entropy*invfmGeV
 end


 
@inline function τ_shear(T,entropy,y::QGPViscosity{N}) where {N}
    viscosity(T,entropy,y)/((T*entropy)*y.Cs)
end


@inline function viscosity(T,entropy,y::ZeroViscosity) 
    zero(promote_type(typeof(T),typeof(entropy)))
 end


 
@inline function τ_shear(T,entropy,y::ZeroViscosity) 
   one(promote_type(typeof(T),typeof(entropy)))
end



function normalization(T,μ,x::Heavy_Quark)
    norm = 0
    for i in x.hadron_list.particle_list
        m = i.Mass
        deg = i.Degeneracy
        q = i.Nc + i.Nac
        n = free_hadron(T,μ,deg,q; m = m)[1]
        norm += q^2*n
    end 
    return norm;    
end




function τ_diffusion(T,x::ZeroDiffusion) 
    return 1.0
end

function τ_diffusion(T,n, x::ZeroDiffusion) 
    return 1.0
end

function τ_diffusion_hadron(T,μ,x::Heavy_Quark,y::ZeroDiffusion) 
    return 1.0
end

function τ_diffusion_hadron(T,n,μ,x::Heavy_Quark,y::ZeroDiffusion) 
    return 1.0
end


function diffusion(T,n,x::ZeroDiffusion) 
    return 0.0
end

function diffusion(T,x::ZeroDiffusion) 
    return 0.0
end


function diffusion_hadron(T,n,μ,x::Heavy_Quark,y::ZeroDiffusion) 
    return 0.0
end

function diffusion_hadron(T,μ,x::Heavy_Quark,y::ZeroDiffusion) 
    return 0.0
end

# shear viscosity and τ_shear of QGPViscosity from a Thermodynamic
@inline function viscosity(T,x::Thermodynamic{N,1,1},y::QGPViscosity{N}) where {N}
    y.ηs*(@inbounds x.pressure_derivative[1])*invfmGeV
end
@inline function τ_shear(T,x::Thermodynamic{N,1,1},y::QGPViscosity{N}) where {N}
    y.ηs*invfmGeV/(T*y.Cs)
end


# every charm-EOS input of the generated charm matrices at (T, α): n and its derivatives from the list, and
# from its species sums weighted by q² n_i the enthalpy per particle h, ⟨m²⟩ and C (Boltzmann: C_i = m K4/K3)
function charm_eos(T, α, list::CharmHadronList)
    th = thermodynamic(T, α, list)
    n = th.pressure
    dndT, dnda = th.pressure_derivative
    R = float(promote_type(typeof(T), typeof(α)))
    # common scale e^{-sh}, sh = max_i (q_i α - z_i), so no weight underflows
    sh = typemin(R)
    for p in list.particle_list
        sh = max(sh, round(p.Nc + p.Nac) * α - p.Mass / T)
    end
    D = zero(R); dD = zero(R); Hn = zero(R); dHn = zero(R); M2 = zero(R); WC = zero(R); dWC = zero(R)
    for p in list.particle_list
        q = round(p.Nc + p.Nac); m = p.Mass; g = p.Degeneracy
        z = m / T
        K1 = besselk1x(z); K2 = besselkx(2, z)          # exponentially scaled; the e^z cancels in ratios
        K3 = K1 + 4 / z * K2
        K4 = K2 + 6 / z * K3
        K5 = K3 + 8 / z * K4
        ex = exp(q * α - z - sh)
        w  = q^2 * g * m^2 * T * K2 * ex                                # q² n_i, up to a common factor
        dw = q^2 * g * ex * m^2 * (m * K1 + 2 * T * K2 + m * K3) / (2 * T)   # ∂_T at fixed α
        hi  = m * K3 / K2
        dhi = -(z / T) * m * ((-(K2 + K4) / 2) * K2 - K3 * (-(K1 + K3) / 2)) / K2^2
        Ci  = m * K4 / K3
        dCi = -(z / T) * m * ((-(K3 + K5) / 2) * K3 - K4 * (-(K2 + K4) / 2)) / K3^2
        D += w; dD += dw
        Hn += w * hi; dHn += dw * hi + w * dhi
        M2 += w * m^2
        WC += w * hi * Ci; dWC += (dw * hi + w * dhi) * Ci + w * hi * dCi
    end
    h = Hn / D; dhdT = (dHn - h * dD) / D
    C = WC / Hn; dCdT = (dWC - C * dHn) / Hn
    return (n = n, dndT = dndT, dnda = dnda, h = h, dhdT = dhdT,
            I31 = n * T * h, dI31dT = dndT * T * h + n * h + n * T * dhdT, dI31da = dnda * T * h,
            aQ = M2 / D + 5 * T * h, bQ = T * h, dbQdT = h + T * dhdT, cQ = C, dcQdT = dCdT)
end

"""
    ColdEdge(; T_edge = 0.156, width = 0.020, floor = 0.01, corona = true)
The cold-edge treatment of the two-way system. Charm decouples below T_edge: D_s → f(T) D_s with f = 1
for T ≥ T_edge, f = floor below T_edge - width, and a C² smoothstep in between. With `corona`,
`two_way_problem` also updates the charm number below T_edge in flux form (the corona update), which
keeps n positive at an under-resolved charm edge.
"""
struct ColdEdge
    T_edge::Float64
    width::Float64
    floor::Float64
    corona::Bool
end
function ColdEdge(; T_edge = 0.156, width = 0.020, floor = 0.01, corona::Bool = true)
    (0 < floor <= 1 && width > 0) || error("ColdEdge: need 0 < floor ≤ 1 and width > 0")
    return ColdEdge(T_edge, width, floor, corona)
end

# the factor f(T) on D_s; edge = nothing switches it off
@inline edge_factor(T, ::Nothing) = one(T)
@inline function edge_factor(T, e::ColdEdge)
    T < e.T_edge || return one(T)
    x = (T - (e.T_edge - e.width)) / e.width
    x <= 0 && return oftype(one(T), e.floor)
    return e.floor + (1 - e.floor) * (x^3 * (10 - 15x + 6x^2))
end

# the weight of the corona update: 0 for T ≥ T_edge, 1 below T_edge - width, the same smoothstep between
@inline function corona_weight(T, e::ColdEdge)
    (e.corona && T < e.T_edge) || return zero(T)
    x = (T - (e.T_edge - e.width)) / e.width
    x <= 0 && return one(T)
    return 1 - x^3 * (10 - 15x + 6x^2)
end



