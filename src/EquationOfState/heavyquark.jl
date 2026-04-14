Base.@kwdef struct Heavy_Quark{T,S} <:EquationOfState 
    mass::T = 1.5
    eos::FluiduMEoS{T} = FluiduMEoS()
    hadron_list::HadronResonaceGas_ccbar{S,T} = HadronResonaceGas_ccbar(readresonancelist(), 30.)
end

function Heavy_Quark(particle_list, ccbar)
    list = HadronResonaceGas_ccbar(readresonancelist(;name_file = particle_list), ccbar)
    return Heavy_Quark(hadron_list = list)
end

@inline @fastmath pressure(T, x::Heavy_Quark) = pressure(T, x.eos)
@inline pressure_derivative(T, v::Val{1}, x::Heavy_Quark) = pressure_derivative(T, v, x.eos)
@inline pressure_derivative(T, v::Val{2}, x::Heavy_Quark) = pressure_derivative(T, v, x.eos)

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

@inline function DsT(m::LinearDiffusion,T; Tc = 0.180)
    if T > Tc
        return m.slope*T + m.offset
    else
        return m.slope*Tc + m.offset
    end
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
        b2 = besselkx(2,m/T)*exp(-m/T)
        b1 = besselk1x(m/T)*exp(-m/T)
        b3 = b1+4/(m/T)*b2  
        b4 = b2 + 6/(m/T)*b3  
        b5 = b3+8/(m/T)*b4
        ex = exp(q*μ)
        tauq += ((2*π*DsT(y,T))/(192*π^3*T^3)*q^2*m^5*ex*(2*b1 - 3*b3 +b5)); 
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



