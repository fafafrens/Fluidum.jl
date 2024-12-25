
struct ThermodynamicPerturbation{T,N,M,G} 
    pressure::T
    pressure_derivative::NTuple{N,T}
    pressure_hessian::NTuple{M,T}
    pressure_third::NTuple{G,T}
end

#just ways how to show thermodynamics. Calling println I will automatically call this Base.show
function Base.show(io::IO, z::ThermodynamicPerturbation{T,N,M,G} ) where {T,N,M,G}
    print(io,"p(x_1,..,x_$N)" , z.pressure, "∇p(x_1,..,x_$N)" ,pressure_derivative,"∇²p(x_1,..,x_$N)" ,pressure_hessian,"∇³p(x_1,..,x_$N)" ,pressure_third  )
end

function Base.show(io::IO, ::MIME"text/plain", z::ThermodynamicPerturbation{T,N,M,G} ) where {T,N,M,G}
    print(io,"p(x_1,..,x_$N)=" , z.pressure,"\n")
    print(io,"∇p(x_1,..,x_$N)=" ,z.pressure_derivative,"\n" )
    print(io,"∇²p(x_1,..,x_$N)=" ,z.pressure_hessian ,"\n" )
    print(io,"∇³p(x_1,..,x_$N)=" ,z.pressure_third ,"\n" )
end

function Base.show(io::IO, z::ThermodynamicPerturbation{T,1,1,1} ) where {T}
    print(io,"p(T)=" , z.pressure, " ∇p(T)=" ,z.pressure_derivative[1]," ∇²p(T)=" ,z.pressure_hessian[1] ,
    " ∇³p(T)=" ,z.pressure_third[1] )
end


function Base.show(io::IO, ::MIME"text/plain", z::ThermodynamicPerturbation{T,1,1,1} ) where {T}
    print(io,"p(T)=" , z.pressure,"\n")
    print(io,"∇p(T)=" ,z.pressure_derivative[1],"\n" )
    print(io,"∇²p(T)=" ,z.pressure_hessian[1] ,"\n" )
    print(io,"∇³p(T)=" ,z.pressure_third[1] ,"\n" )
end
