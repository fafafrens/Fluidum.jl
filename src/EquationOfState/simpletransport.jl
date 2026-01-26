struct SimpleShearViscosity{T} <:ShearViscosity
    ηs::T
    Cs::T
end 

struct ZeroViscosity<:ShearViscosity

end

@inline function viscosity(T,μ,x::Thermodynamic{N,2,3},y::SimpleShearViscosity{N}) where {N}
   entropy=x.pressure_derivative[1]
    y.ηs*entropy*invfmGeV
end
@inline function viscosity(T,μ,x::Thermodynamic{N,2,3},y::ZeroViscosity) where {N}
    zero(promote_type(typeof(T),typeof(μ)))
end

@inline function viscosity(T,x::Thermodynamic{N,1,1},y::ZeroViscosity)  where {N}
    zero(promote_type(typeof(T)))
end

@inline function viscosity(T,x::ThermodynamicPerturbation{N,1,1,1},y::ZeroViscosity)  where {N}
    zero(promote_type(typeof(T),N))
end


@inline function viscosity_derivative(T,x::ThermodynamicPerturbation{N,1,1,1},y::ZeroViscosity)  where {N}
    zero(promote_type(typeof(T),N))
end

@inline function viscosity(T,x::Thermodynamic{N,1,1},y::SimpleShearViscosity{N}) where {N}
    entropy= @inbounds x.pressure_derivative[1]
     y.ηs*entropy*invfmGeV
 end

 @inline function viscosity(T,x::ThermodynamicPerturbation{N,1,1,1},y::SimpleShearViscosity{N}) where {N}
    entropy= @inbounds x.pressure_derivative[1]
     y.ηs*entropy*invfmGeV
 end

 @inline function viscosity_derivative(T,x::ThermodynamicPerturbation{N,1,1,1},y::SimpleShearViscosity{N}) where {N}
    dentropy= @inbounds x.pressure_hessian[1]
     y.ηs*dentropy*invfmGeV
 end


@inline function τ_shear(T,μ,x::Thermodynamic{N,2,3},y::SimpleShearViscosity{N}) where {N}
    viscosity(T,μ,x,y)/((T*x.pressure_derivative[1]+μ*x.pressure_derivative[2])*y.Cs)
end

@inline function τ_shear(T,x::Thermodynamic{N,1,1},y::SimpleShearViscosity{N}) where {N}
    #entropy=x.pressure_derivative[1]
     y.ηs*invfmGeV/(T*y.Cs)
    
    #viscosity(T,x,y)/((T*x.pressure_derivative[1])*y.Cs)
end

@inline function τ_shear(T,x::ThermodynamicPerturbation{N,1,1,1},y::SimpleShearViscosity{N}) where {N}
    #entropy=x.pressure_derivative[1]
     y.ηs*invfmGeV/(T*y.Cs)
    
    #viscosity(T,x,y)/((T*x.pressure_derivative[1])*y.Cs)
end

@inline function τ_shear_dervative(T,x::ThermodynamicPerturbation{N,1,1,1},y::SimpleShearViscosity{N}) where {N}
    #entropy=x.pressure_derivative[1]
     -y.ηs*invfmGeV/(T^2*y.Cs)
    
    #viscosity(T,x,y)/((T*x.pressure_derivative[1])*y.Cs)
end
@inline function τ_shear_dervative(T,x::ThermodynamicPerturbation{N,1,1,1},y::ZeroViscosity) where {N}
    #entropy=x.pressure_derivative[1]
    zero(promote_type(typeof(T),N))
    #viscosity(T,x,y)/((T*x.pressure_derivative[1])*y.Cs)
end
@inline function τ_shear(T,μ,x::Thermodynamic{N,2,3},y::ZeroViscosity)  where {N}
    one(promote_type(typeof(T),typeof(μ)))
end

@inline function τ_shear(T,x::Thermodynamic{N,1,1},y::ZeroViscosity)  where {N}
    one(promote_type(typeof(T)))
end

@inline function τ_shear(T,μ,x::ThermodynamicPerturbation{N,1,1,1},y::ZeroViscosity) where {N}
    #entropy=x.pressure_derivative[1]
    one(promote_type(typeof(T),typeof(μ)))
    
    #viscosity(T,x,y)/((T*x.pressure_derivative[1])*y.Cs)
end

@inline function τ_shear(T,x::ThermodynamicPerturbation{N,1,1,1},y::ZeroViscosity) where {N}
    #entropy=x.pressure_derivative[1]
    one(promote_type(typeof(T)))
end


struct SimpleDiffusionCoefficient{T}<:Diffusion
    κ::T
    Cκ::T
end

struct ZeroDiffusion{T}<:Diffusion
    mass::T
end

ZeroDiffusion(; mass::T = 1.5) where {T<:Real} = ZeroDiffusion{T}(mass)
ZeroDiffusion(mass::T) where {T<:Real} = ZeroDiffusion{T}(mass)



@inline function diffusion(T,x::Thermodynamic{N,1,1},y::SimpleDiffusionCoefficient{N}) where{N}
    zero(typeof(T))
end

@inline function diffusion(T,x::ThermodynamicPerturbation{N,1,1,1},y::SimpleDiffusionCoefficient{N}) where{N}
    zero(typeof(T))
end

@inline function diffusion(T,μ,x::Thermodynamic{N,2,3},y::SimpleDiffusionCoefficient{N}) where{N}
    fmGeV*fmGeV*T^4/μ^2*y.κ/( 2*π)     #y.κ*(μ^2)*( T*x.pressure_derivative[2]   )^2/(T*x.pressure_derivative[1]+μ*x.pressure_derivative[2])^2*invfmGeV^2
end

@inline function τ_diffusion(T,μ,x::Thermodynamic{N,2,3},y::SimpleDiffusionCoefficient{N}) where{N}
    1/y.Cκ*diffusion(T,μ,x,y)/x.pressure_hessian[3]*(μ+T*(T*x.pressure_hessian[1]+μ*x.pressure_hessian[2])/(T*x.pressure_hessian[2]+μ*x.pressure_hessian[3]))*1/T^2
end

@inline function τ_diffusion(T,x::Thermodynamic{N,1,1},y::SimpleDiffusionCoefficient{N}) where{N}
    one(typeof(T))
end
@inline function τ_diffusion(T,x::ThermodynamicPerturbation{N,1,1,1},y::SimpleDiffusionCoefficient{N}) where{N}
    one(typeof(T))
end


@inline function diffusion(T,μ,x::Thermodynamic{N,2,3},y::ZeroDiffusion)  where {N}
    zero(promote_type(typeof(T),typeof(μ)))
end

@inline function diffusion(T,x::Thermodynamic{N,1,1},y::ZeroDiffusion)  where {N}
    zero(promote_type(typeof(T)))
end

@inline function diffusion(T,x::ThermodynamicPerturbation{N,1,1,1},y::ZeroDiffusion)  where {N}
    zero(promote_type(typeof(T)))
end

@inline function τ_diffusion(T,μ,x::Thermodynamic{N,2,3},y::ZeroDiffusion)  where {N}
    one(promote_type(typeof(T),typeof(μ)))
end

@inline function τ_diffusion(T,x::Thermodynamic{N,1,1},y::ZeroDiffusion)  where {N}
    one(promote_type(typeof(T)))
end

@inline function τ_diffusion(T,x::ThermodynamicPerturbation{N,1,1,1},y::ZeroDiffusion)  where {N}
    one(promote_type(typeof(T)))
end

struct SimpleBulkViscosity{T} <:BulkViscosity
    ζs::T
    Cζ::T
end


struct ZeroBulkViscosity<:BulkViscosity
end



#@inline @fastmath function speed_of_sound_squared(T,x::EquationOfState)
#    pressure_derivative(T,Val(1),x)/(T*pressure_derivative(T,Val(2),x)) 
#end

@inline function bulk_viscosity(T,x::Thermodynamic{N,1,1},y::SimpleBulkViscosity{N}) where {N}
   y.ζs/(1+((T-0.175)/0.024)^2)*invfmGeV*x.pressure_derivative[1]
end
@inline function τ_bulk(T,x::Thermodynamic{N,1,1},y::SimpleBulkViscosity{N}) where {N}
    cs2= x.pressure_derivative[1]/(T*x.pressure_hessian[1])
    bulk_viscosity(T,x,y)/(T*x.pressure_derivative[1]*y.Cζ)*1/(1/3-cs2)^2+0.1
end



@inline function bulk_viscosity(T,x::ThermodynamicPerturbation{N,1,1,1},y::SimpleBulkViscosity{N}) where {N}
    y.ζs/(1+((T-0.175)/0.024)^2)*invfmGeV*x.pressure_derivative[1]
 end

 @inline function τ_bulk(T,x::ThermodynamicPerturbation{N,1,1,1},y::SimpleBulkViscosity{N}) where {N}
     cs2= x.pressure_derivative[1]/(T*x.pressure_hessian[1])
     bulk_viscosity(T,x,y)/(T*x.pressure_derivative[1]*y.Cζ)*1/(1/3-cs2)^2+0.1
 end


@inline function bulk_viscosity(T,α,x::Thermodynamic{N,2,3},y::SimpleBulkViscosity{N}) where {N}
    y.ζs/(1+((T-0.175)/0.024)^2)*invfmGeV*x.pressure_derivative[1]
end

@inline function τ_bulk(T, α,x::Thermodynamic{N,2,3},y::SimpleBulkViscosity{N}) where {N}
    cs2=0.2
    bulk_viscosity(T,α,x,y)/(T*x.pressure_derivative[1]*y.Cζ)*1/(1/3-cs2)^2
end

@inline function bulk_viscosity(T,μ,x::Thermodynamic{N,2,3},y::ZeroBulkViscosity)  where {N}
    zero(promote_type(typeof(T),typeof(μ)))
end

@inline function bulk_viscosity(T,x::Thermodynamic{N,1,1},y::ZeroBulkViscosity)  where {N}
    zero(promote_type(typeof(T)))
end

@inline function bulk_viscosity(T,x::ThermodynamicPerturbation{N,1,1,1},y::ZeroBulkViscosity)  where {N}
    zero(promote_type(typeof(T),N))
end

@inline function bulk_viscosity_derivative(T,x::ThermodynamicPerturbation{N,1,1,1},y::ZeroBulkViscosity)  where {N}
    zero(promote_type(typeof(T),N))
end



@inline function τ_bulk(T,μ,x::Thermodynamic{N,2,3},y::ZeroBulkViscosity)  where {N}
    one(promote_type(typeof(T),typeof(μ)))
end

@inline function τ_bulk(T,x::Thermodynamic{N,1,1},y::ZeroBulkViscosity)  where {N}
    one(promote_type(typeof(T)))
end

@inline function τ_bulk(T,x::ThermodynamicPerturbation{N,1,1,1},y::ZeroBulkViscosity)  where {N}
    one(promote_type(typeof(T),N))
end

@inline function τ_bulk_derivative(T,x::ThermodynamicPerturbation{N,1,1,1},y::ZeroBulkViscosity)  where {N}
    zero(promote_type(typeof(T),N))
end