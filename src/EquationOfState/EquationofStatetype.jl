

abstract type EquationOfState end

abstract type CompositeEquationOfState <:EquationOfState end

iscomposite(::T) where {T<:EquationOfState}=false
iscomposite(::T) where {T<:CompositeEquationOfState} =true
    
                        

struct SumEquationOfState{A <:EquationOfState ,B <:EquationOfState} <:CompositeEquationOfState 
    equation_of_state_1::A
    equation_of_state_2::B
end 
function Base.show(io::IO, z::SumEquationOfState) 
    print(io,"(",z.equation_of_state_1,")","+","(",z.equation_of_state_2,")"  )
end


struct DifferenceEquationOfState{A <:EquationOfState ,B <:EquationOfState} <:CompositeEquationOfState 
    equation_of_state_1::A
    equation_of_state_2::B
end
function Base.show(io::IO, z::DifferenceEquationOfState) 
    print(io,"(",z.equation_of_state_1,")","-","(",z.equation_of_state_2,")"  )
end

struct ProductEquationOfState{A <:EquationOfState ,B <:EquationOfState} <:CompositeEquationOfState 
    equation_of_state_1::A
    equation_of_state_2::B
end
function Base.show(io::IO, z::ProductEquationOfState) 
    print(io,"(",z.equation_of_state_1,")","*","(",z.equation_of_state_2,")"  )
end

struct DivisionEquationOfState{A <:EquationOfState ,B <:EquationOfState} <:CompositeEquationOfState 
    equation_of_state_1::A
    equation_of_state_2::B
end
function Base.show(io::IO, z::DivisionEquationOfState) 
    print(io,"(",z.equation_of_state_1,")","/","(",z.equation_of_state_2 ,")" )
end

struct OppositeEquationOfState{B <:EquationOfState} <:CompositeEquationOfState 
    equation_of_state::B
end 
function Base.show(io::IO, z::OppositeEquationOfState) 
    print(io,"-","(",z.equation_of_state,")"  )
end

struct InverseEquationOfState{B <:EquationOfState} <:CompositeEquationOfState 
    equation_of_state::B
end 
function Base.show(io::IO, z::InverseEquationOfState) 
    print(io,"1/","(",z.equation_of_state,")"  )
end

struct NumberSumEquationOfState{A <:Number,B <:EquationOfState} <:CompositeEquationOfState 
    num::A
    equation_of_state_2::B
end
function Base.show(io::IO, z::NumberSumEquationOfState) 
    print(io,"(",z.num,")","+","(",z.equation_of_state_2,")"  )
end


struct NumberDifferenceEquationOfState{A <:Number ,B <:EquationOfState} <:CompositeEquationOfState 
    num::A
    equation_of_state_2::B
end
function Base.show(io::IO, z::NumberDifferenceEquationOfState) 
    print(io,"(",z.num,")","-","(",z.equation_of_state_2,")"  )
end

struct NumberProductEquationOfState{A <:Number ,B <:EquationOfState} <:CompositeEquationOfState 
    num::A
    equation_of_state_2::B
end
function Base.show(io::IO, z::NumberProductEquationOfState) 
    print(io,"(",z.num,")","*","(",z.equation_of_state_2,")"  )
end

struct NumberDivisionEquationOfState{A <:Number ,B <:EquationOfState} <:CompositeEquationOfState 
    num::A
    equation_of_state_2::B
end
function Base.show(io::IO, z::NumberDivisionEquationOfState) 
    print(io,"(",z.num,")","/","(",z.equation_of_state_2,")"  )
end

Base.:+(a::A,b::B) where {A<:EquationOfState,B<:EquationOfState}=SumEquationOfState{A,B}(a,b)
Base.:-(a::A,b::B) where {A<:EquationOfState,B<:EquationOfState}=DifferenceEquationOfState{A,B}(a,b)
Base.:*(a::A,b::B) where {A<:EquationOfState,B<:EquationOfState}=ProductEquationOfState{A,B}(a,b)
Base.:/(a::A,b::B) where {A<:EquationOfState,B<:EquationOfState}=DivisionEquationOfState{A,B}(a,b)

Base.:+(a::A,b::B) where {A<:EquationOfState,B<:Number}=NumberSumEquationOfState{B,A}(b,a)
Base.:-(a::A,b::B) where {A<:EquationOfState,B<:Number}=NumberSumEquationOfState{B,A}(-b,a)
Base.:*(a::A,b::B) where {A<:EquationOfState,B<:Number}=NumberProductEquationOfState{B,A}(b,a)
function Base.:/(a::A,b::B) where {A<:EquationOfState,B<:Number}
    binv=inv(b)
    NumberProductEquationOfState{promote_type(B,typeof(binv)),A}(binv,a)
end
Base.:+(a::A,b::B) where {A<:Number,B<:EquationOfState}=NumberSumEquationOfState{A,B}(a,b)
Base.:-(a::A,b::B) where {A<:Number,B<:EquationOfState}=NumberSumEquationOfState{A,OppositeEquationOfState{B}}(a,OppositeEquationOfState{B}(b))
Base.:*(a::A,b::B) where {A<:Number,B<:EquationOfState}=NumberProductEquationOfState{A,B}(a,b)
#implement the inv(a) on the equation of state
Base.:/(a::A,b::B) where {A<:Number,B<:EquationOfState}=NumberProductEquationOfState{A,B}(a,inv(b))
Base.:-(b::B) where {B<:EquationOfState}=OppositeEquationOfState{B}(b)




 thermodynamic(T,x::SumEquationOfState)= thermodynamic(T,x.equation_of_state_1)+thermodynamic(T,x.equation_of_state_2)
 thermodynamic(T,x::DifferenceEquationOfState)= thermodynamic(T,x.equation_of_state_1)-thermodynamic(T,x.equation_of_state_2)
 thermodynamic(T,x::ProductEquationOfState)= thermodynamic(T,x.equation_of_state_1)*thermodynamic(T,x.equation_of_state_2)
 thermodynamic(T,x::DivisionEquationOfState)= thermodynamic(T,x.equation_of_state_1)/thermodynamic(T,x.equation_of_state_2)
 
 thermodynamic(T,x::NumberSumEquationOfState)= x.num+thermodynamic(T,x.equation_of_state_2)
 thermodynamic(T,x::NumberDifferenceEquationOfState)= x.num-thermodynamic(T,x.equation_of_state_2)
 thermodynamic(T,x::NumberProductEquationOfState)= x.num*thermodynamic(T,x.equation_of_state_2)
 thermodynamic(T,x::NumberDivisionEquationOfState)= x.num/thermodynamic(T,x.equation_of_state_2)
 thermodynamic(T,x::OppositeEquationOfState)= -thermodynamic(T,x.equation_of_state)
 thermodynamic(T,x::InverseEquationOfState)= 1/thermodynamic(T,x.equation_of_state)


 pressure(T,wal::N) where {N<:EquationOfState}= thermodynamic(T,wal).pressure

pressure_derivative(T,::Val{1},wal::N) where {N<:EquationOfState}= thermodynamic(T,wal).pressure_derivative[1]
pressure_derivative(T,::Val{2},wal::N) where {N<:EquationOfState}= thermodynamic(T,wal).pressure_hessian[1]
pressure_derivative(T,::Val{3},wal::N) where {N<:EquationOfState}= thermodynamic_perturbation(T,wal).pressure_third[1]
entropy(T,wal::N) where {N<:EquationOfState}= thermodynamic(T,wal).pressure_derivative[1]





 thermodynamic(T,mu,x::SumEquationOfState)= thermodynamic(T,mu,x.equation_of_state_1)+thermodynamic(T,mu,x.equation_of_state_2)
 thermodynamic(T,mu,x::DifferenceEquationOfState)= thermodynamic(T,mu,x.equation_of_state_1)-thermodynamic(T,mu,x.equation_of_state_2)
 thermodynamic(T,mu,x::ProductEquationOfState)= thermodynamic(T,mu,x.equation_of_state_1)*thermodynamic(T,mu,x.equation_of_state_2)
 thermodynamic(T,mu,x::DivisionEquationOfState)= thermodynamic(T,mu,x.equation_of_state_1)/thermodynamic(T,mu,x.equation_of_state_2)
 
 thermodynamic(T,mu,x::NumberSumEquationOfState)= x.num+thermodynamic(T,mu,x.equation_of_state_2)
 thermodynamic(T,mu,x::NumberDifferenceEquationOfState)= x.num-thermodynamic(T,mu,x.equation_of_state_2)
 thermodynamic(T,mu,x::NumberProductEquationOfState)= x.num*thermodynamic(T,mu,x.equation_of_state_2)
 thermodynamic(T,mu,x::NumberDivisionEquationOfState)= x.num/thermodynamic(T,mu,x.equation_of_state_2)
 thermodynamic(T,mu,x::OppositeEquationOfState)= -thermodynamic(T,mu,x.equation_of_state)



 pressure(T,mu,wal::N) where {N<:EquationOfState}= thermodynamic(T,mu,wal).pressure

 pressure_derivative(T,mu,::Val{1},::Val{0},wal::N) where {N<:EquationOfState}= thermodynamic(T,mu,wal).pressure_derivative[1]
 pressure_derivative(T,mu,::Val{0},::Val{1},wal::N) where {N<:EquationOfState}= thermodynamic(T,mu,wal).pressure_derivative[2]
 pressure_derivative(T,mu,::Val{2},::Val{0},wal::N) where {N<:EquationOfState}= thermodynamic(T,mu,wal).pressure_hessian[1]
 pressure_derivative(T,mu,::Val{1},::Val{1},wal::N) where {N<:EquationOfState}= thermodynamic(T,mu,wal).pressure_hessian[2]
 pressure_derivative(T,mu,::Val{0},::Val{2},wal::N) where {N<:EquationOfState}= thermodynamic(T,mu,wal).pressure_hessian[3]
    
function energy_density(T,mu,wal::N) where {N<:EquationOfState}
    x=thermodynamic(T,mu,wal)
     -x.pressure+T*x.pressure_derivative[1]+mu*x.pressure_derivative[2]
end

 struct OneDPicewiseEquationOfState{A<:EquationOfState,B<:AbstractInterval ,C<:EquationOfState,D<:AbstractInterval} <:CompositeEquationOfState
    equation_of_state_1::A
    range_1::B
    equation_of_state_2::C
    range_2::D

    function OneDPicewiseEquationOfState(eos1::A,int1::B,eos2::C,int2::D) where {A,B,C,D}
        intersect=IntervalSets.intersect(int1, int2)
        l,r=endpoints(intersect)

        if isempty(intersect)
            #if they are not 
            return new{A,B,C,D}(eos1,int1,eos2,int2)
        end 

        if isempty(intersect)==false && l==r
            ##here we shrink one of the interval
            return new{A,B,C,D}(eos1,int1,eos2,int2)
        end 

        if isempty(intersect)==false && l!=r

        return throw(ArgumentError("The interval are not disjoint i really do not know what to do."))

        end 

        
    end
 end

 function Base.show(io::IO, z::OneDPicewiseEquationOfState)
    print(io,z.equation_of_state_1," for T ⊆ ",z.range_1," and ",z.equation_of_state_2," for T ⊆ ",z.range_2 )
end


function ⊕(x::Tuple{A,B},y::Tuple{C,D}) where {A<:EquationOfState,B<:AbstractInterval,C<:EquationOfState,D<:AbstractInterval}
    OneDPicewiseEquationOfState(x...,y...)
end

function ⊕(x::Tuple{A,B,C},y::Tuple{D,E,F}) where {A<:EquationOfState,B<:AbstractInterval,C<:AbstractInterval,D<:EquationOfState,E<:AbstractInterval,F<:AbstractInterval}
    TwoDPicewiseEquationOfState(x...,y...)
end


function ⊕(x::Tuple{A,Tuple{B,C}},y::Tuple{D,Tuple{E,F}}) where {A<:EquationOfState,B<:AbstractInterval,C<:AbstractInterval,D<:EquationOfState,E<:AbstractInterval,F<:AbstractInterval}
    TwoDPicewiseEquationOfState(x...,y...)
end



struct TwoDPicewiseEquationOfState{A<:EquationOfState,B<:Tuple{AbstractInterval,AbstractInterval},C<:EquationOfState,D<:Tuple{AbstractInterval,AbstractInterval}} <:CompositeEquationOfState
    equation_of_state_1::A
    range_1::B
    equation_of_state_2::C
    range_2::D

    function TwoDPicewiseEquationOfState(eos1::A,int1::B,eos2::C,int2::D) where {A,B,C,D}#improve border check
        intersect=IntervalSets.intersect.(int1, int2)
        boundary_intersect=endpoints.(intersect)
        
        if !all(isempty.(intersect))&&(isequal(boundary_intersect[2][1],boundary_intersect[2][2])||
            isequal(boundary_intersect[1][1],boundary_intersect[1][2])
            )
            
            return new{A,B,C,D}(eos1,int1,eos2,int2)
        
        end

        if all(isempty.(intersect))
            #if they are not 
            return new{A,B,C,D}(eos1,int1,eos2,int2)
        end 

         

        if !all(isempty.(intersect)) && !isequal(boundary_intersect[1][1],boundary_intersect[1][2])&&!isequal(boundary_intersect[2][1],boundary_intersect[2][2])

        return throw(ArgumentError("The interval are not disjoint i really do not know what to do."))

        end 

        
    end

 end


 function Base.show(io::IO, z::TwoDPicewiseEquationOfState)
    print(io,z.equation_of_state_1," for T,μ ⊆ ",z.range_1[1]," ⊕ ",z.range_1[2],", ",z.equation_of_state_2," for T,μ ⊆ ",z.range_2[1]," ⊕ ",z.range_2[2])
end


function Base.show(io::IO, ::MIME"text/plain", z::TwoDPicewiseEquationOfState)
    print(io,z.equation_of_state_1," for T,μ ⊆ ",z.range_1[1]," ⊕ ",z.range_1[2],"\n")
    print(io,z.equation_of_state_2," for T,μ ⊆ ",z.range_2[1]," ⊕ ",z.range_2[2])
end

 function TwoDPicewiseEquationOfState(eos1::A,int1::B, int2::C ,eos2::D,int3::E,int4::F) where{A<:EquationOfState,B<:AbstractInterval,C<:AbstractInterval,D<:EquationOfState,E<:AbstractInterval,F<:AbstractInterval}

    TwoDPicewiseEquationOfState(eos1,(int1,int2),eos2,(int3,int4)) 
 end



 function thermodynamic(T::S,x::OneDPicewiseEquationOfState{A,B,C,D})::Thermodynamic{S,1,1} where {S<:Number,A<:EquationOfState,B<:AbstractInterval,C<:EquationOfState,D<:AbstractInterval} 
    if T ∈  x.range_1
        return thermodynamic(T,x.equation_of_state_1)
    end 
    
    if T ∈  x.range_2
        return thermodynamic(T,x.equation_of_state_2)
    end 

    return Thermodynamic{S,1,1}(0,(0,),(0,))
end


function thermodynamic(T::S,mu::G,x::TwoDPicewiseEquationOfState{A,B,C,D})::Thermodynamic{promote_type(S,G),2,3} where {S<:Number,G<:Number,A<:EquationOfState,B,C<:EquationOfState,D} 
    
    if all( (T, mu) .∈  x.range_1 ) 
        return thermodynamic(T,mu,x.equation_of_state_1)
    end 
    
    if all( (T, mu) .∈  x.range_2 ) 
        return thermodynamic(T,mu,x.equation_of_state_2)
    end 

    return Thermodynamic{promote_type(S,G),2,3}(0,(0,0),(0,0,0))
end

    



