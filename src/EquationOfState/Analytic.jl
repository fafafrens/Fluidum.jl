
struct OneDAnalyticEquationOfState{D,D1,D2} <:EquationOfState
    f::D
    f_1::D1
    f_2::D2
    
end

function Base.show(io::IO, z::OneDAnalyticEquationOfState)
    print(io,z.f,"(T)")
end


struct TwoDAnalyticEquationOfState{F,F10,F01,F20,F11,F02} <:EquationOfState
    f::F
    f_10::F10
    f_01::F01
    f_20::F20
    f_11::F11
    f_02::F02
end 

function Base.show(io::IO, z::TwoDAnalyticEquationOfState)
    print(io,z.f,"(T,μ)")
end

Analytic(f::A,f_1::B,f_2::C) where {A,B,C}=OneDAnalyticEquationOfState(f,f_1,f_2)

Analytic(f::A,f_10::B,f_01::C,f_20::D,f_11::E,f_02::F) where {A,B,C,D,E,F}=TwoDAnalyticEquationOfState(f,f_10,f_01,f_20,f_11,f_02)


struct Gluing{A} <:EquationOfState
    f::A
end

function Base.show(io::IO, z::Gluing{A}) where A
    print(io,z.f)
end

#function Gluing(f::T) where {T}
#    Gluing{T}(f)
#end

thermodynamic(T::A,x::OneDAnalyticEquationOfState) where{A<:Number} = Thermodynamic{A,1,1}(x.f(T),(x.f_1(T),),(x.f_2(T),))

thermodynamic(T::A,mu::B,x::TwoDAnalyticEquationOfState) where{A<:Number,B<:Number} = Thermodynamic{promote_type(A,B),2,3}(x.f(T,mu),(x.f_10(T,mu),x.f_01(T,mu)),(x.f_20(T,mu),x.f_11(T,mu),x.f_02(T,mu)))

thermodynamic(T::A,x::Gluing) where{A<:Number} = Thermodynamic{A,1,1}(x.f(T),(zero(A),),(zero(A),))

function thermodynamic(T::A,mu::B,x::Gluing) where{A<:Number,B<:Number}
    N=promote_type(A,B)
    Thermodynamic{N,2,3}(x.f(T,mu),(zero(N),zero(N)),(zero(N),zero(N),zero(N)))
end
