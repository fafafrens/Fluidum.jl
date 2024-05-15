
abstract type AbstractBoundedFunction{left,right} end


@inline leftbounds(::AbstractBoundedFunction{left,right}) where {left,right} =left 
@inline rightbounds(::AbstractBoundedFunction{left,right} ) where {left,right} =right 
@inline isinbounds(x,::AbstractBoundedFunction{left,right}) where {left,right}=  all(left .<= x .<=right )

struct BoundedFunction{X,left,right} <: AbstractBoundedFunction{left,right}
    func::X
    
    function BoundedFunction(f::X,left::NTuple{N,S},right::NTuple{N,T}) where {X,N,S,T}
        new{typeof(f),left,right}(f)
    end
end 

BoundedFunction(f::X,left::T,right::S) where {X,T <:Number,S<:Number}=BoundedFunction(f, ( left,),(right, ))

@inline leftbounds(::BoundedFunction{X,left,right}) where {X,left,right} =left 
@inline rightbounds(::BoundedFunction{X,left,right}) where {X,left,right} =right 
@inline isinbounds(x,::BoundedFunction{X,left,right}) where {X,left,right}=  all(left .<= x .<=right )

@inline (f::BoundedFunction{X,left,right})(x::NTuple{N,T}) where {X,left,right,N,T} =f.func(x)
@inline (f::BoundedFunction{X,left,right})(x::SVector{N,T}) where {X,left,right,N,T} =f.func(ntuple(i->x[i],Val{N}()))


struct PieceWiseFunction{X,Y,left,right} <: AbstractBoundedFunction{left,right}
    fun1::X
    fun2::Y
  
end

function PieceWiseFunction(f::X,g::Y) where {X<: AbstractBoundedFunction,Y <: AbstractBoundedFunction }
        
    left=min.(leftbounds(f),leftbounds(g))
    right=max.(rightbounds(f),rightbounds(g))
    PieceWiseFunction{typeof(f),typeof(g),left,right}(f,g)
end

function PieceWiseFunction(f::X) where {X<: AbstractBoundedFunction}
        
    left=leftbounds(f)
    right=rightbounds(f)
    zerofunction(x)=zero(f(left))
    g=BoundedFunction(zerofunction,left.*Inf,right.*(-Inf))

    PieceWiseFunction{typeof(f),typeof(g),left,right}(f,g)
end

function PieceWiseFunction(f::X,g::Y,h::Z) where {X<: AbstractBoundedFunction,Y<: AbstractBoundedFunction,Z<: AbstractBoundedFunction}
        
    rest=PieceWiseFunction(g,h)
    PieceWiseFunction(f,rest)
end


function PieceWiseFunction(f...)
    
    PieceWiseFunction(first(f),PieceWiseFunction(Base.tail(f)...))
end

PieceWiseFunction(f::X,left::NTuple{N,S},right::NTuple{N,T}) where {X,N,S,T}=PieceWiseFunction(BoundedFunction(f,left,right))

PieceWiseFunction(f::X,left::S,right::T) where {X,S<:Number,T<:Number}=PieceWiseFunction(BoundedFunction(f,(left,),(right,)))

@inline leftbounds(::PieceWiseFunction{X,Y,left,right}) where {X,Y,left,right} =left 
@inline rightbounds(::PieceWiseFunction{X,Y,left,right}) where {X,Y,left,right} =right 
@inline isinbounds(x,::PieceWiseFunction{X,Y,left,right}) where {X,Y,left,right}=  all(left .<= x .<=right )



(fun::PieceWiseFunction{X,Y,left,right})(x...) where {X,Y,left,right}=fun(x)


function (fun::PieceWiseFunction{X,Y,left,right})(x::NTuple{N,T}) where {X,Y,left,right,N,T}
    if isinbounds(x,fun.fun1) 
        return fun.fun1(x)
    elseif isinbounds(x,fun.fun2) 
        return fun.fun2(x)
    else
        return  zero(first(promote(
                fun.fun1(leftbounds(fun.fun1)), fun.fun2(leftbounds(fun.fun2))
            )))
    end 
    
end

struct SplineInterp{tuple,A,left,right} <:AbstractBoundedFunction{left,right}
a::A
end

function SplineInterp(fun,tuple)
    left=leftbounds(fun)
    right=rightbounds(fun)
    ranges=range.(leftbounds(fun),rightbounds(fun),tuple)
    iter =Iterators.product(ranges...)

    A = [fun(i) for i in iter ]

    itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
    
    sitp = scale(itp, ranges...)

    SplineInterp{tuple,typeof(sitp),left,right}(sitp)
end 

function (f::SplineInterp{tuple,A,left,right})(x::NTuple{N,T}) where {tuple,A,left,right,N,T}
    f.a(x...)
end

function (f::SplineInterp{tuple,A,left,right})(x::SVector{N,T}) where {tuple,A,left,right,N,T}
    f.a(x...)
end

function (f::SplineInterp{tuple,A,left,right})(x...) where {tuple,A,left,right}
    f.a(x...)
end


jacobian(f::SplineInterp{tuple,A,left,right},x::NTuple{N,T}) where {tuple,A,left,right,N,T} =reduce(hcat,Interpolations.gradient(f.a,x...))'

jacobian(f::SplineInterp{tuple,A,left,right},x::SVector{N,T}) where {tuple,A,left,right,N,T}  =reduce(hcat,Interpolations.gradient(f.a,x...))'

jacobian(f::SplineInterp{tuple,A,left,right},x...) where {tuple,A,left,right}  =reduce(hcat,Interpolations.gradient(f.a,x...))'