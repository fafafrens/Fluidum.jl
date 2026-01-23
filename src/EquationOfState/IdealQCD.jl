struct IdealQCD{T} <:EquationOfState
    a1::T
    a2::T
    a3::T
end

function IdealQCD(Nc,Nf)  
    a1=8*pi^2/15 *(Nc^2 -1 +7/4* Nc*Nf) 
    a2= 2.0*Nc *Nf / 27
    a3= 2*Nc *Nf / (81*pi^2)
    IdealQCD(a1,a2,a3) 
end

IdealQCD() =IdealQCD(3,2)

function IdealQCDT(Nc,Nf)
    a1=8*pi^2/15 *(Nc^2 -1 +7/4* Nc*Nf) 
    a2= 0.
    a3= 0.
    IdealQCD(a1,a2,a3) 
    
end
#StronglyCoupled()+IdealQCD()



@inline  function pressure(T,x::IdealQCD) 
    1/(4*3*2) *x.a1 *T^4 *fmGeV3
end

@inline  function pressure(T,μ,x::IdealQCD) 
    fmGeV3*(1/(4*3*2) *x.a1 *T^4 +1/4*x.a2 * T^2*μ^2 +1/(4*3*2) *x.a3* μ^4)
end

@inline  function pressure_derivative(T,μ,::Val{1},::Val{0},x::IdealQCD ) 
    (1/(3*2) *x.a1 *T^3 +2/4*x.a2 * T*μ^2 )*fmGeV3
end

@inline  function pressure_derivative(T,::Val{1},x::IdealQCD ) 
    1/(3*2) *x.a1 *T^3 *fmGeV3
end

@inline @fastmath function pressure_derivative(T,μ,::Val{2},::Val{0},x::IdealQCD ) 
    
   ( 1/2 *x.a1 *T^2 +2/4*x.a2 *μ^2 )*fmGeV3
end

@inline @fastmath function pressure_derivative(T,μ,::Val{3},::Val{0},x::IdealQCD ) 
    
    ( x.a1 *T )*fmGeV3
end


@inline @fastmath function pressure_derivative(T,::Val{2},x::IdealQCD ) 
    
    1/2 *x.a1 *T^2 *fmGeV3
end

@inline @fastmath function pressure_derivative(T,::Val{3},x::IdealQCD ) 
    
     x.a1 *T *fmGeV3
end

@inline @fastmath function pressure_derivative(T,μ,::Val{0},::Val{1},x::IdealQCD ) 
    
    (2/4*x.a2 * T^2*μ +1/(3*2) *x.a3* μ^3 )*fmGeV3
end

@inline @fastmath function pressure_derivative(T,μ,::Val{0},::Val{2},x::IdealQCD ) 
    
    (2/4*x.a2 * T^2 +1/(2) *x.a3* μ^2)*fmGeV3
end

@inline @fastmath function pressure_derivative(T,μ,::Val{1},::Val{1},x::IdealQCD ) 
    (x.a2 * T*μ )*fmGeV3
end


thermodynamic(T::N,x::IdealQCD{S}) where {N,S} = Thermodynamic{promote_type(N,S),1,1}(pressure(T,x),(pressure_derivative(T,Val{1}(),x),),(pressure_derivative(T,Val{2}(),x),) ) 

@inline @fastmath function thermodynamic_perturbation(T,x::IdealQCD{<:Number})
    

    
    ThermodynamicPerturbation(pressure(T,x),(pressure_derivative(T,Val{1}(),x),),(pressure_derivative(T,Val{2}(),x),) ,(pressure_derivative(T,Val{3}(),x),))

end




function thermodynamic(T::N,α::S,x::IdealQCD) where {N,S}
    Thermodynamic{promote_type(N,S),2,3}(pressure(T,α,x),
    (pressure_derivative(T,α,Val{1}(),Val{0}(),x),pressure_derivative(T,α,Val{0}(),Val{1}(),x)), 
    (pressure_derivative(T,α,Val{2}(),Val{0}(),x),pressure_derivative(T,α,Val{1}(),Val{1}(),x),pressure_derivative(T,α,Val{0}(),Val{2}(),x) ))
end


