


abstract type TransportCoefficient end

struct FluidProperties{T,S}
    equation_of_state::T
    transport_coefficient::S
end



struct NoTransportProperties<:TransportCoefficient

end 


struct Conformal{T} <:TransportCoefficient
    ηs::T
    κ::T
    Cs::T
    Cκ::T
end

#StronglyCoupled()=Conformal(1/4/pi,4-2*log(2),pi/2)
#StronglyCoupled()=Conformal(1/4/pi,,4-2*log(2),pi/2)
#StronglyCoupled()=Conformal(1/4/pi,1/4/pi,4-2*log(2),pi/4)

struct FullViscous{T} <:TransportCoefficient
    ηs::T
    κ::T
    ζs::T
    Cs::T
    Cκ::T
    Cζ::T
end



@inline viscosity(T,μ,x::EquationOfState,y::Conformal)=y.ηs*entropy_density(T,μ,x)*invfmGeV
@inline viscosity(T,μ,x::EquationOfState,y::FullViscous)=y.ηs*entropy_density(T,μ,x)*invfmGeV
@inline viscosity(T,μ,x::EquationOfState,y::NoTransportProperties)=0

@inline viscosity(T,x::EquationOfState,y::Conformal)=y.ηs*entropy_density(T,x)*invfmGeV
@inline viscosity(T,x::EquationOfState,y::FullViscous)=y.ηs*entropy_density(T,x)*invfmGeV
@inline viscosity(T,x::EquationOfState,y::NoTransportProperties)=0

@inline τ_shear(T,μ,x::EquationOfState,y::Conformal)=viscosity(T,μ,x,y)/((T*entropy_density(T,μ,x)+μ*charge_density(T,μ,x))*y.Cs) 
@inline τ_shear(T,μ,x::EquationOfState,y::FullViscous)= viscosity(T,μ,x,y)/((T*entropy_density(T,μ,x)+μ*charge_density(T,μ,x))*y.Cs) 
@inline τ_shear(T,μ,x::EquationOfState,y::NoTransportProperties)=1

@inline τ_shear(T,x::EquationOfState,y::Conformal)= viscosity(T,x,y)/(T*entropy_density(T,x)*y.Cs)
@inline τ_shear(T,x::EquationOfState,y::FullViscous)= viscosity(T,x,y)/(T*entropy_density(T,x)*y.Cs)
@inline τ_shear(T,x::EquationOfState,y::NoTransportProperties)= 1

@inline bulk_viscosity(T,μ,x::EquationOfState,y::Conformal)=0
@inline bulk_viscosity(T,μ,x::EquationOfState,y::FullViscous)=0#y.ζs/(1+((T-0.175)/0.024)^2)*fmGeV*entropy_density(T,μ,x)
@inline bulk_viscosity(T,μ,x::EquationOfState,y::NoTransportProperties)=0

@inline bulk_viscosity(T,x::EquationOfState,y::Conformal)=0
@inline bulk_viscosity(T,x::EquationOfState,y::FullViscous)=y.ζs/(1+((T-0.175)/0.024)^2)*fmGeV*entropy_density(T,x)
@inline bulk_viscosity(T,x::EquationOfState,y::NoTransportProperties)=0

@inline τ_bulk(T,μ,x::EquationOfState,y::Conformal)=1
@inline τ_bulk(T,μ,x::EquationOfState,y::FullViscous)=1#bulk_viscosity(T,μ,x,y)/((T*entropy_density(T,μ,x)+μ*charge_density(T,μ,x))*y.Cζ)*1/(1/3-speed_of_sound_squared(T,μ,x))
@inline τ_bulk(T,μ,x::EquationOfState,y::NoTransportProperties)=1

@inline τ_bulk(T,x::EquationOfState,y::Conformal)=1
@inline τ_bulk(T,x::EquationOfState,y::FullViscous)=bulk_viscosity(T,x,y)/(T*entropy_density(T,x)*y.Cζ)*1/(1/3-speed_of_sound_squared(T,x))
@inline τ_bulk(T,x::EquationOfState,y::NoTransportProperties)=1

#@inline diffusion(T,μ,x::EquationOfState,y::Conformal)= y.κ*( charge_density(T,μ,x))^2/(T*entropy_density(T,μ,x)+μ*charge_density(T,μ,x))^2*invfmGeV
#@inline diffusion(T,μ,x::EquationOfState,y::FullViscous)= y.κ*( charge_density(T,μ,x))^2/(T*entropy_density(T,μ,x)+μ*charge_density(T,μ,x))^2*invfmGeV
#@inline diffusion(T,μ,x::EquationOfState,y::NoTransportProperties)= 1
#@inline diffusion(T,μ,x::FluidProperties)=diffusion(T,μ,x.equation_of_state,x.transport_coefficient)

@inline diffusion(T,μ,x::EquationOfState,y::Conformal)= y.κ*(T^2+T*μ+μ^2)*( T*charge_density(T,μ,x))^2/(T*entropy_density(T,μ,x)+μ*charge_density(T,μ,x))^2*invfmGeV^2
@inline diffusion(T,μ,x::EquationOfState,y::FullViscous)= y.κ*(μ^2)*( T*charge_density(T,μ,x))^2/(T*entropy_density(T,μ,x)+μ*charge_density(T,μ,x))^2*invfmGeV^2
@inline diffusion(T,μ,x::EquationOfState,y::NoTransportProperties)= 1


#@inline function τ_diffusion(T,μ,x::EquationOfState,y::Conformal)
 #   1y.τκ*diffusion(T,μ,x,y)/(pressure_derivative(T,μ,0,2,x))
    #/(pressure_derivative(T,μ,0,2,x))*
    #(T*(T*pressure_derivative(T,μ,2,0,x)+μ*pressure_derivative(T,μ,1,1,x))/(T*pressure_derivative(T,μ,1,1,x)+μ*pressure_derivative(T,μ,0,2,x))+μ)
#end
#@inline function τ_diffusion(T,μ,x::EquationOfState,y::FullViscous)
 #   y.τκ*diffusion(T,μ,x,y)/(pressure_derivative(T,μ,0,2,x))
 #   /(pressure_derivative(T,μ,0,2,x))*
 #   (T*(T*pressure_derivative(T,μ,2,0,x)+μ*pressure_derivative(T,μ,1,1,x))/(T*pressure_derivative(T,μ,1,1,x)+μ*pressure_derivative(T,μ,0,2,x))+μ)
#end

#@inline function τ_diffusion(T,μ,x::EquationOfState,y::NoTransportProperties)
#    1y.τκ*diffusion(T,μ,x,y)/(pressure_derivative(T,μ,0,2,x))
#    #/(pressure_derivative(T,μ,0,2,x))*
#    #(T*(T*pressure_derivative(T,μ,2,0,x)+μ*pressure_derivative(T,μ,1,1,x))/(T*pressure_derivative(T,μ,1,1,x)+μ*pressure_derivative(T,μ,0,2,x))+μ)
#end

@inline function τ_diffusion(T,μ,x::EquationOfState,y::FullViscous)
    1/y.Cκ*diffusion(T,μ,x,y)/pressure_derivative(T,μ,Val(0),Val(2),x)*(μ+T*(T*pressure_derivative(T,μ,Val(2),Val(0),x)+μ*pressure_derivative(T,μ,Val(1),Val(1),x))/(T*pressure_derivative(T,μ,Val(1),Val(1),x)+μ*pressure_derivative(T,μ,Val(0),Val(2),x)))*1/T^2
end
    
@inline function τ_diffusion(T,μ,x::EquationOfState,y::NoTransportProperties)
    1#1/y.Cκ*diffusion(T,μ,x,y)/pressure_derivative(T,μ,Val(0),Val(2),x)*(μ+T*(T*pressure_derivative(T,μ,Val(2),Val(0),x)+μ*pressure_derivative(T,μ,Val(1),Val(1),x))/(T*pressure_derivative(T,μ,Val(1),Val(1),x)+μ*pressure_derivative(T,μ,Val(0),Val(2),x)))*1/T^2
end








#lamda for everything execpt T and mu 0
#@inline @fastmath function lambda(T,μ,x::IdealQCD)
#    (pressure_derivative(T,μ,Val(1),Val(0),x)^2*pressure_derivative(T,μ,Val(0),Val(2),x)+pressure_derivative(T,μ,Val(0),Val(1),x)*(pressure_derivative(T,μ,Val(2),Val(0),x)*pressure_derivative(T,μ,Val(1),Val(0),x)-2*pressure_derivative(T,μ,Val(1),Val(0),x)*pressure_derivative(T,μ,Val(1),Val(1),x)))/((pressure_derivative(T,μ,Val(2),Val(0),x)*pressure_derivative(T,μ,Val(0),Val(2),x)-pressure_derivative(T,μ,Val(1),Val(1),x)^2)*(T*pressure_derivative(T,μ,Val(1),Val(0),x)+μ*pressure_derivative(T,μ,Val(0),Val(1),x)))
#end

#@inline lambda(T,μ,x::EquationOfState)=lambda(T,μ,x)
#@inline lambda(T,μ,x::FluidProperties)= lambda(T,μ,x.equation_of_state)