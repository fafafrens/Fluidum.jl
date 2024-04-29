

abstract type TransportCoefficient end

struct FluidProperties{T,S} where {T<:EquationOfState,S<:TransportCoefficient}
    equation_of_state::T
    transport_coefficient::S
end

import Base.+
+(x::EquationOfState,y::TransportCoefficient)=FluidProperties(x,y)
+(y::TransportCoefficient,x::EquationOfState)=FluidProperties(x,y)

@inline  pressure(T,μ,x::FluidProperties) = pressure(T,μ,x.equation_of_state) 

@inline  pressure(T,x::FluidProperties) = pressure(T,x.equation_of_state) 


pressure_derivative(T,μ,::Val{2},::Val{0},x::FluidProperties )=pressure_derivative(T,μ,Val(2),Val(0),x.equation_of_state)

pressure_derivative(T,μ,::Val{1},::Val{0},x::FluidProperties )=pressure_derivative(T,μ,Val(1),Val(0),x.equation_of_state)


pressure_derivative(T,::Val{1},x::FluidProperties )=pressure_derivative(T,Val(1),x.equation_of_state)

pressure_derivative(T,::Val{2},x::FluidProperties )=pressure_derivative(T,Val(2),x.equation_of_state)

pressure_derivative(T,μ,::Val{0},::Val{1},x::FluidProperties )=pressure_derivative(T,μ,Val(0),Val(1),x.equation_of_state)


pressure_derivative(T,μ,::Val{0},::Val{2},x::FluidProperties )=pressure_derivative(T,μ,Val(0),Val(2),x.equation_of_state)

pressure_derivative(T,μ,::Val{1},::Val{1},x::FluidProperties )=pressure_derivative(T,μ,Val(1),Val(1),x.equation_of_state)




@inline @fastmath function energy_density(T,x::EquationOfState)
    pressure_derivative(T,Val(1),x)*T - pressure(T,x)  
end



@inline @fastmath function energy_density(T,μ,x::EquationOfState)
   -pressure(T,μ,x)+T*pressure_derivative(T,μ,Val(1),Val(0),x)+μ*pressure_derivative(T,μ,Val(0),Val(1),x)
end

@inline energy_density(T,μ,x::FluidProperties)= energy_density(T,μ,x.equation_of_state)

# T derivative of energy density
@inline @fastmath function energy_density_derivative(T,μ,::Val{1},::Val{0},x::EquationOfState)
    T*pressure_derivative(T,μ,Val(2),Val(0),x) + μ*pressure_derivative(T,μ,Val(1),Val(1),x)
end

# T derivative of energy density (specific heat) why mu below but not in def?
@inline @fastmath function energy_density_derivative(T,::Val{1},x::EquationOfState)
    T*pressure_derivative(T,μ,Val(2),Val(0),x)
end

# mu derivative of energy density
@inline @fastmath function energy_density_derivative(T,μ,::Val{0},::Val{1},x::EquationOfState)
    T*pressure_derivative(T,μ,Val(1),Val(1),x) + pressure_derivative(T,μ,Val(0),Val(2),x)*μ
end

@inline @fastmath function speed_of_sound_squared(T,x::EquationOfState)
    pressure_derivative(T,Val(1),x)/(T*pressure_derivative(T,Val(2),x)) 
end

#physical speed of sound (at fixed s/n)
@inline @fastmath function speed_of_sound_squared(T,μ,x::EquationOfState)
    (pressure_derivative(T,μ,Val(0),Val(1),x)^2*pressure_derivative(T,μ,Val(2),Val(0),x)-2*pressure_derivative(T,μ,Val(1),Val(0),x)*pressure_derivative(T,μ,Val(0),Val(1),x)*pressure_derivative(T,μ,Val(1),Val(1),x)+pressure_derivative(T,μ,Val(1),Val(0),x)^2*pressure_derivative(T,μ,Val(0),Val(2),x))/((T*pressure_derivative(T,μ,Val(1),Val(0),x)+μ*pressure_derivative(T,μ,Val(0),Val(1),x))*(pressure_derivative(T,μ,Val(2),Val(0),x)*pressure_derivative(T,μ,Val(0),Val(2),x)-pressure_derivative(T,μ,Val(1),Val(1),x)^2)) 
end

#@inline pressure_derivative(T,μ,n,m,x::EquationOfState) =pressure_derivative(T,μ,Val{n}(),Val{m}(),x ) 
#@inline pressure_derivative(T,μ,n,m,x::FluidProperties) =pressure_derivative(T,μ,Val{n}(),Val{m}(),x.equation_of_state ) 
     
@inline function pressure_derivative(T,μ,n::Int,m::Int,x::EquationOfState)
    if n==1 
        if m==1 
            return pressure_derivative(T,μ,Val(1),Val(1),x )
        elseif m==2
            return pressure_derivative(T,μ,Val(1),Val(2),x )
        else 
            throw(ArgumentError("The derivative are implement up to the second order. If you really want it 
        you could implemnt pressure_derivative(T,μ,Val(n),Val(m),x::EquationOfState)  for the sepcific type of equation of state or 
        use AD on this function "))
        end 
    elseif n==2
        if m==1 
            return pressure_derivative(T,μ,Val(2),Val(1),x )
        elseif m==2
            return pressure_derivative(T,μ,Val(2),Val(2),x )
        else 
            throw(ArgumentError("The derivative are implement up to the second order. If you really want it 
        you could implemnt pressure_derivative(T,μ,Val(n),Val(3),x::EquationOfState)  for the sepcific type of equation of state or 
        use AD on this function "))
        end
    else 
        throw(ArgumentError("The derivative are implement up to the second order. If you really want it 
        you could implemnt pressure_derivative(T,μ,Val(3),Val(m),x::EquationOfState)  for the sepcific type of equation of state or 
        use AD on this function "))
    end 

    
end
@inline pressure_derivative(T,μ,n::Int,m::Int,x::FluidProperties) =pressure_derivative(T,μ,n,m,x.equation_of_state ) 

@inline function pressure_derivative(T,n::Int,x::EquationOfState)
    if n==1 
        return pressure_derivative(T,Val(1),x )
    elseif n==2
        return pressure_derivative(T,Val(2),x )
    else 
        throw(ArgumentError("The derivative are implement up to the second order. If you really want it 
        you could implemnt pressure_derivative(T,Val(3),x::EquationOfState) and the successive for the sepcific type of equation of state or 
        use AD on this function "))
    end 
   
end
@inline pressure_derivative(T,n::Int,x::FluidProperties) =pressure_derivative(T,n,x.equation_of_state ) 



@inline entropy_density(T,μ,x::EquationOfState)= pressure_derivative(T,μ,Val(1),Val(0),x)
@inline entropy_density(T,μ,x::FluidProperties)= entropy_density(T,μ,x.equation_of_state)
@inline entropy_density(T,x::EquationOfState)= pressure_derivative(T,Val(1),x)
@inline entropy_density(T,x::FluidProperties)= entropy_density(T,x.equation_of_state)


@inline charge_density(T,μ,x::EquationOfState)= pressure_derivative(T,μ,Val(0),Val(1),x)
@inline charge_density(T,μ,x::FluidProperties)= charge_density(T,μ,x.equation_of_state)


@inline viscosity(T,μ,x::FluidProperties)= viscosity(T,μ,x.equation_of_state,x.transport_coefficient)

@inline viscosity(T,x::FluidProperties)= viscosity(T,x.equation_of_state,x.transport_coefficient)


@inline τ_shear(T,μ,x::FluidProperties)=τ_shear(T,μ,x.equation_of_state,x.transport_coefficient)

@inline τ_shear(T,x::FluidProperties)=τ_shear(T,x.equation_of_state,x.transport_coefficient)


@inline bulk_viscosity(T,μ,x::FluidProperties)=bulk_viscosity(T,μ,x.equation_of_state,x.transport_coefficient)

@inline bulk_viscosity(T,x::FluidProperties)=bulk_viscosity(T,x.equation_of_state,x.transport_coefficient)

@inline τ_bulk(T,μ,x::FluidProperties)=τ_bulk(T,μ,x.equation_of_state,x.transport_coefficient)

@inline τ_bulk(T,x::FluidProperties)=τ_bulk(T,x.equation_of_state,x.transport_coefficient)

@inline diffusion(T,μ,x::FluidProperties)=diffusion(T,μ,x.equation_of_state,x.transport_coefficient)

@inline τ_diffusion(T,μ,x::FluidProperties)=τ_diffusion(T,μ,x.equation_of_state,x.transport_coefficient)


thermodyanmics(T,mu,wal::FluidProperties) = thermodyanmics(T,mu,wal.equation_of_state)
