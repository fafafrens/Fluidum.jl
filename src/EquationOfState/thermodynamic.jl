
struct Thermodynamic{T,N,M} 
    pressure::T
    pressure_derivative::NTuple{N,T}
    pressure_hessian::NTuple{M,T}
end

#just ways how to show thermodynamics. Calling println I will automatically call this Base.show
function Base.show(io::IO, z::Thermodynamic{T,N,M} ) where {T,N,M}
    print(io,"p(x_1,..,x_$N)" , z.pressure, "∇p(x_1,..,x_$N)" ,pressure_derivative,"∇²p(x_1,..,x_$N)" ,pressure_hessian  )
end

function Base.show(io::IO, ::MIME"text/plain", z::Thermodynamic{T,N,M} ) where {T,N,M}
    print(io,"p(x_1,..,x_$N)=" , z.pressure,"\n")
    print(io,"∇p(x_1,..,x_$N)=" ,z.pressure_derivative,"\n" )
    print(io,"∇²p(x_1,..,x_$N)=" ,z.pressure_hessian ,"\n" )
end

function Base.show(io::IO, z::Thermodynamic{T,1,1} ) where {T}
    print(io,"p(T)=" , z.pressure, " ∇p(T)=" ,z.pressure_derivative[1]," ∇²p(T)=" ,z.pressure_hessian[1]  )
end


function Base.show(io::IO, ::MIME"text/plain", z::Thermodynamic{T,1,1}) where {T}
    print(io,"p(T)=" , z.pressure,"\n")
    print(io,"∇p(T)=" ,z.pressure_derivative[1],"\n" )
    print(io,"∇²p(T)=" ,z.pressure_hessian[1] ,"\n" )
end

function Base.show(io::IO, z::Thermodynamic{T,2,3} ) where {T}
    print(io,"p(T,μ)=" , z.pressure, " ∇p(T,μ)=" ,z.pressure_derivative," ∇²p(T,μ)=" ,z.pressure_hessian  )
end

function Base.show(io::IO, ::MIME"text/plain", z::Thermodynamic{T,2,3}) where{T}
    print(io,"p(T,μ)=" , z.pressure,"\n")
    print(io,"∇p(T,μ)=" ,z.pressure_derivative,"\n" )
    print(io,"∇²p(T,μ)=" ,z.pressure_hessian ,"\n" )
end

#way to define the sum between two stuctures of thermodynamics type
@inline function Base.:+(a::Thermodynamic{T,N,M},b::Thermodynamic{T,N,M}) where {T,N,M} 
    Thermodynamic{T,N,M}(
    a.pressure+b.pressure,
    ntuple(i->a.pressure_derivative[i]+b.pressure_derivative[i],Val(N)),
    ntuple(i->a.pressure_hessian[i]+b.pressure_hessian[i],Val(M)))
end

@inline function Base.:+(a::Thermodynamic{T,N,M},b::Thermodynamic{S,N,M}) where {T,N,M,S} 
        Thermodynamic{promote_type(T,S),N,M}(
        a.pressure+b.pressure,
        ntuple(i->a.pressure_derivative[i]+b.pressure_derivative[i],Val(N)),
        ntuple(i->a.pressure_hessian[i]+b.pressure_hessian[i],Val(M)))
end


@inline function Base.:-(a::Thermodynamic{T,N,M},b::Thermodynamic{T,N,M}) where {T,N,M} 
    Thermodynamic{T,N,M}(
    a.pressure-b.pressure,
    ntuple(i->a.pressure_derivative[i]-b.pressure_derivative[i],Val(N)),
    ntuple(i->a.pressure_hessian[i]-b.pressure_hessian[i],Val(M)))
end

@inline function Base.:-(a::Thermodynamic{T,N,M}) where {T,N,M} 
    Thermodynamic{T,N,M}(
    -a.pressure,
    ntuple(i->-a.pressure_derivative[i],Val(N)),
    ntuple(i->-a.pressure_hessian[i],Val(M)))
end


@inline function Base.:-(a::Thermodynamic{T,N,M},b::Thermodynamic{S,N,M}) where {T,N,M,S} 
    Thermodynamic{promote_type(T,S),N,M}(
    a.pressure-b.pressure,
    ntuple(i->a.pressure_derivative[i]-b.pressure_derivative[i],Val(N)),
    ntuple(i->a.pressure_hessian[i]-b.pressure_hessian[i],Val(M)))
end



@inline function Base.:+(a::S,b::Thermodynamic{S,N,M}) where {S,N,M} 
    Thermodynamic{S,N,M}(
    a+b.pressure,
    ntuple(i->b.pressure_derivative[i],Val(N)),
    ntuple(i->b.pressure_hessian[i],Val(M)))
end

@inline function Base.:+(a::T,b::Thermodynamic{S,N,M}) where {T,S,N,M} 
    Thermodynamic{promote_type(T,S),N,M}(
    a+b.pressure,
    ntuple(i->b.pressure_derivative[i],Val(N)),
    ntuple(i->b.pressure_hessian[i],Val(M)))
end





 @inline function Base.:-(a::S,b::Thermodynamic{S,N,M}) where {S,N,M} 
    Thermodynamic{S,N,M}(
    a-b.pressure,
    ntuple(i->-b.pressure_derivative[i],Val(N)),
    ntuple(i->-b.pressure_hessian[i],Val(M)))
end

@inline function Base.:-(a::T,b::Thermodynamic{S,N,M}) where {T,S,N,M} 
    Thermodynamic{promote_type(T,S),N,M}(
    a-b.pressure,
    ntuple(i->-b.pressure_derivative[i],Val(N)),
    ntuple(i->-b.pressure_hessian[i],Val(M)))
end

@inline Base.:-(a::Thermodynamic{S,N,M},b::Number) where {S,N,M} = a+(-b)

 @inline function Base.:*(a::S,b::Thermodynamic{S,N,M}) where {S,N,M} 
    Thermodynamic{S,N,M}(
    a*b.pressure,
    ntuple(i->a*b.pressure_derivative[i],Val(N)),
    ntuple(i->a*b.pressure_hessian[i],Val(M)))
end

@inline function Base.:*(a::T,b::Thermodynamic{S,N,M}) where {T,S,N,M} 
    Thermodynamic{promote_type(T,S),N,M}(
    a*b.pressure,
    ntuple(i->a*b.pressure_derivative[i],Val(N)),
    ntuple(i->a*b.pressure_hessian[i],Val(M)))
end

@inline Base.:+(a::Thermodynamic{S,N,M},b::Number) where {S,N,M} = b*a

@inline Base.:/(a::Number,b::Thermodynamic{S,N,M}) where {S,N,M} = a*inv(b)

@inline Base.:/(a::Thermodynamic{S,N,M},b::Number) where {S,N,M} = a*inv(b)
    


    

 #lindx(i,j,n)=(-(-1 + i)*(i - 2* n))÷2+j
 
 #axes(StaticArrays.SHermitianCompact{2}(SVector{3}(1,2,3)))
 

@inline function Base.:*(a::Thermodynamic{T,1,1},b::Thermodynamic{T,1,1}) where {T} 
    Thermodynamic{T,1,1}(a.pressure*b.pressure,
    (a.pressure_derivative[1]*b.pressure+b.pressure_derivative[1]*a.pressure,),
    (a.pressure_hessian[1]*b.pressure+2*a.pressure_derivative[1]*b.pressure_derivative[1]+b.pressure_hessian[1]*a.pressure,))
end

@inline function Base.:*(a::Thermodynamic{S,1,1},b::Thermodynamic{T,1,1}) where {S,T} 
    Thermodynamic{promote_type(T,S),1,1}(a.pressure*b.pressure,
    (a.pressure_derivative[1]*b.pressure+b.pressure_derivative[1]*a.pressure,),
    (a.pressure_hessian[1]*b.pressure+2*a.pressure_derivative[1]*b.pressure_derivative[1]+b.pressure_hessian[1]*a.pressure,))
end

@inline function Base.inv(a::Thermodynamic{T,1,1}) where {T}
    Thermodynamic{T,1,1}(1/a.pressure,(-a.pressure_derivative[1]/a.pressure^2,),
    (2*a.pressure_derivative[1]^2/a.pressure^3-a.pressure_hessian[1]^2/a.pressure^2 ,))

end 


@inline function Base.:*(a::Thermodynamic{T,2,3},b::Thermodynamic{T,2,3}) where {T}
    Thermodynamic{T,2,3}(a.pressure*b.pressure,
        (a.pressure_derivative[1]*b.pressure+b.pressure_derivative[1]*a.pressure,
            a.pressure_derivative[2]*b.pressure+b.pressure_derivative[2]*a.pressure ),
        ( a.pressure_hessian[1]*b.pressure+2*a.pressure_derivative[1]*b.pressure_derivative[1]+b.pressure_hessian[1]*a.pressure,
            a.pressure_derivative[1]*b.pressure_derivative[2]+a.pressure_derivative[2]*b.pressure_derivative[1]
                +b.pressure_hessian[2]*a.pressure +a.pressure_hessian[2]*b.pressure,
            a.pressure_hessian[3]*b.pressure+2*a.pressure_derivative[2]*b.pressure_derivative[2]+b.pressure_hessian[3]*a.pressure)
    )
    
    
end

@inline function Base.:*(a::Thermodynamic{S,2,3},b::Thermodynamic{T,2,3}) where {S,T}
    Thermodynamic{promote_type(T,S),2,3}(a.pressure*b.pressure,
        (a.pressure_derivative[1]*b.pressure+b.pressure_derivative[1]*a.pressure,
            a.pressure_derivative[2]*b.pressure+b.pressure_derivative[2]*a.pressure ),
        ( a.pressure_hessian[1]*b.pressure+2*a.pressure_derivative[1]*b.pressure_derivative[1]+b.pressure_hessian[1]*a.pressure,
            a.pressure_derivative[1]*b.pressure_derivative[2]+a.pressure_derivative[2]*b.pressure_derivative[1]
                +b.pressure_hessian[2]*a.pressure +a.pressure_hessian[2]*b.pressure,
            a.pressure_hessian[3]*b.pressure+2*a.pressure_derivative[2]*b.pressure_derivative[2]+b.pressure_hessian[3]*a.pressure)
    )
    
    
end

@inline function Base.inv(a::Thermodynamic{T,2,3}) where {T}

    Thermodynamic{T,2,3}(1/a.pressure,
        (
        -a.pressure_derivative[1]/a.pressure^2 
        ,-a.pressure_derivative[2]/a.pressure^2 
        )
        ,(2*a.pressure_derivative[1]^2/a.pressure^3-a.pressure_hessian[1]^2/a.pressure^2 , 
            2*a.pressure_derivative[1]*a.pressure_derivative[2]/a.pressure^3-a.pressure_hessian[2]/a.pressure^2
            ,2*a.pressure_derivative[2]^2/a.pressure^3-a.pressure_hessian[3]^2/a.pressure^2)
        )


end 

@inline Base.:/(a::Thermodynamic{T,N,M},b::Thermodynamic{T,N,M}) where {T,N,M}=  a * inv(b)



