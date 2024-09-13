

# this part is fot a given array (Nfield,.....) build up an the information of the index structure 
struct IndexStructure{dimension_tot,dimension_space,N_field,TupleDim}
    left_ghost::Vector{Vector{CartesianIndex{dimension_tot}}}
    right_ghost::Vector{Vector{CartesianIndex{dimension_tot}}}
    interior::Vector{CartesianIndex{dimension_space}}
    unit_inidices::NTuple{dimension_space,CartesianIndex{dimension_tot}}  
    unit_inidices_space::NTuple{dimension_space,CartesianIndex{dimension_space}}  
end

e_i(i,::Val{N}) where {N} = CartesianIndex(ntuple(j-> (j==i ? 1 : 0),Val{N}()  ))

function end_index(i,I::CartesianIndex{N},size) where {N}
    CartesianIndex(ntuple(j-> (j==i ? size : I[j]),Val{N}()  ))
end

function first_index(i,I::CartesianIndex{N}) where {N}
    CartesianIndex(ntuple(j-> (j==i ? 1 : I[j]),Val{N}()  ))
end


"""
    IndexStructure(array)

Create an index structure form an array with ghost cell
"""
function IndexStructure(array)
    sizes=size(array)
    @assert all(sizes[2:end].> 2) "Too few points in the discretization "
    fullindices=CartesianIndices(array)
    left_ghost=[filter(i-> i[j]==1,fullindices) for j in 2:length(sizes) ]
    right_ghost=[filter(i-> i[j]==sizes[j],fullindices) for j in 2:length(sizes) ]
    unit_vector=ntuple(i->e_i(i+1,Val{length(sizes)}()),length(sizes[2:end]))
    unit_vector_space=ntuple(i->e_i(i,Val{length(sizes[2:end])}()),length(sizes[2:end]))
    interior= vec(collect(CartesianIndices(range.(Ref(2),sizes[2:end]))))
    IndexStructure{length(sizes),length(sizes[2:end]),
    sizes[1],sizes[2:end]}(
        left_ghost,right_ghost,interior
        ,unit_vector,unit_vector_space)
end


"""
    IndexStructure(N_field::Int,sizes::NTuple{N,Int}) where {N}

Create an index structure from the size of the field and the space-diemsension. The sizes are the full sizes with ghost cell.
"""
function IndexStructure(N_field::Int,sizes::NTuple{N,Int}) where {N}
    tot_size=(N_field,sizes...)
    tot_dimension=length(tot_size)
    space_dimension=length(sizes)
    tot_ranges=range.(1,tot_size)
    space_ranges=range.(2,sizes .-1)
    
    fullindices=CartesianIndices(tot_ranges)
    
    left_ghost=[filter(i-> i[j]==1,fullindices) for j in range(2,tot_dimension) ]
    right_ghost=[filter(i-> i[j]==tot_size[j],fullindices) for j in range(2,tot_dimension) ]

    unit_vector=ntuple(i->e_i(i+1,Val{tot_dimension}()),space_dimension)
    unit_vector_space=ntuple(i->e_i(i,Val{space_dimension}()),space_dimension)
    interior= vec(collect(CartesianIndices(space_ranges)))
    
    IndexStructure{tot_dimension,space_dimension,
    N_field,sizes}(
        left_ghost,right_ghost,interior
        ,unit_vector,unit_vector_space)


end

#this part handle the generation of the field information 

abstract type DField end

struct NDField{S<:Symbol,N} <:DField
    left_boundary::NTuple{N,S} 
    right_boundary::NTuple{N,S} 
    name::S
end 

#"""
#    NDField(left::NTuple{N,S},right::NTuple{N,S},name::S) where {S<:Symbol}
#
#Create a field with name and left and right boundary information 
#"""
#function NDField(left::NTuple{N,S},right::NTuple{N,S},name::S) where {N,S<:Symbol}
#    NDField(left,right,name)
#end

dimension(::NDField{S,N}) where {S,N} =N
name(x::DField)=string(x.name)
name_symbol(x::DField)=x.name


struct Fields{N_field,N_dimension,S}
    field::NTuple{N_field,NDField{S,N_dimension}}  
end

 

field_dimension(::Fields{N,S,T}) where {N,S,T}= N

names(x::Fields{N,S,T}) where {N,S,T}=name_symbol.(x.field)

get_index(express,x::Fields{N,S,T}) where {N,S,T}=findfirst(x->x==express,names(x))

dimension(::Fields{N,S,T}) where {N,S,T}= S

Fields(x::NDField{S,N_dimension}) where {S,N_dimension}=
Fields{1,N_dimension,S}((x,))

Fields(x::NDField{S,N_dimension},y::NDField{S,N_dimension}) where {S,N_dimension}=
Fields{2,N_dimension,S}((x,y))


Fields(x::Fields{Nf,N,S}) where {Nf,N,S}=x
Fields(x::Fields{Nf,N,S},y::NDField{S,N}) where {Nf,N,S}=Fields(x,Fields(y))
Fields(y::NDField{S,N},x::Fields{Nf,N,S},) where {Nf,N,S}=Fields(Fields(y),x)

Fields(x::Fields{Nf,N,S},y::Fields{Nf2,N,S}) where {Nf,N,S,Nf2}=Fields{Nf+Nf2,N,S}((x.field...,y.field...)) 
#Fields(x,y...)=Fields(Fields(Fields(x),Fields(y[1])),Fields(Base.tail(y)...))
"""
    Fields(y...)

Collect the field all toghther. 
"""
Fields(y...)=Fields(Fields(Base.first(y)),Fields(Base.tail(y)...))





function left_parity(x::NDField{S,N},i::Int) where {S,N}
    #@assert isless(i,N) "Try to access the $i field with number of field $N "
    if x.left_boundary[i]==:odd
        return -1
    elseif x.left_boundary[i]==:even 
        return 1
    elseif x.left_boundary[i]==:ghost
        return 1
    else 
        return 0
    end 
end

function right_parity(x::NDField{S,N},i::Int) where {S,N}
    #@assert isless(i,N) "Try to access the $i field with number of field $N "
    if x.right_boundary[i]==:odd
        return -1
    elseif x.right_boundary[i]==:even 
        return 1
    elseif x.right_boundary[i]==:ghost
        return 1
    else 
        return 0
    end  
end


function periodicity(x::NDField{S,N},i::Int) where {S,N}
    #@assert isless(i,N) "Try to access the $i field with number of field $N "
    if x.right_boundary[i]==x.left_boundary[i]==:periodic
        return 1
    elseif x.right_boundary[i]==x.left_boundary[i]==:antiperiodic
        return -1
    else 
        return 0
    end  

end 


function get_left_parity(Φ::Fields{N_field,N_dimension,S}) where {N_field,N_dimension,S}
   
   
    ntuple(j->SVector{N_field,Int}(ntuple(i->left_parity(Φ.field[i],j),Val{N_field}())),Val{N_dimension}() )
 

 
 end
 
 
 
 function get_right_parity(Φ::Fields{N_field,N_dimension,S}) where {N_field,N_dimension,S}
     
     ntuple(j->SVector{N_field,Int}(ntuple(i->right_parity(Φ.field[i],j),Val{N_field}())),Val{N_dimension}() )
    
 end
 
 function get_periodicity(Φ::Fields{N_field,N_dimension,S}) where {N_field,N_dimension,S}
     ntuple(j->SVector{N_field,Int}(ntuple(i->periodicity(Φ.field[i],j),Val{N_field}())),Val{N_dimension}() )

 end


# this structure hold a chace vector and matrix of the size of the number of field 

struct FieldsCache{T,N_field,N_dim,N_field2}
    vector::Vector{NTuple{N_dim,MVector{N_field,T}}}
    matrix::Vector{NTuple{N_dim,MMatrix{N_field,N_field,T,N_field2}}}
end




function FieldsCache{T,N_field,N_dim,N_field2}(N_copies) where {T,N_field,N_dim,N_field2}
    
    FieldsCache{T,N_field,N_dim,N_field2}( 

[ ntuple(j-> MVector{N_field,T}(undef), Val(N_dim) ) for i in Base.OneTo(N_copies)] 
,
[ ntuple(j-> MMatrix{N_field,N_field,T,N_field2}(undef) , Val(N_dim) ) for i in Base.OneTo(N_copies) ]
)
end

function FieldsCache(type::Type{T},x::Fields{N_field,N_dim,S},N_copies) where {T,N_field,N_dim,S}
    FieldsCache{T,N_field,N_dim,N_field*N_field}(N_copies)
end

function convert_cache(::Type{T},x::FieldsCache{S,N_field,N_dim,N_field2}) where {T,S,N_field,N_dim,N_field2}
    FieldsCache{T,N_field,N_dim,N_field2}(length(x.vector))
end



struct DiscreteFileds{T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,S,N_field2}
    fields::Fields{N_field,space_dimension,S}
    discretization::CartesianDiscretization{space_dimension,Sizes,Lengths,DXS,M}     #Ninterval{space_dimension,Sizes,DXS,M}
    index_structure::IndexStructure{total_dimensions,space_dimension,N_field,Sizes_ghosted}
    left_parity::NTuple{space_dimension,SVector{N_field,Int64}}
    right_parity::NTuple{space_dimension,SVector{N_field,Int64}}
    periodicity::NTuple{space_dimension,SVector{N_field,Int64}}
    cache::FieldsCache{T,N_field,space_dimension,N_field2}
end 


"""
    DiscreteFileds(fields::Fields,discretization::CartesianDiscretization,type::Type{T},N_copies::Int=4) 

Collect the information of the field and the space-doiscretization in one data type. 
"""
function DiscreteFileds(fields::Fields{N_field,space_dimension,S}
,discretization::CartesianDiscretization{space_dimension,Sizes,Lengths,DXS,M},type::Type{T},N_copies::Int=6) where {N_field,space_dimension,S,Sizes,Lengths,DXS,M,T}
if N_copies<6 
    N_copies=6
end 

Sizes_ghosted=Sizes.+2
idx=IndexStructure(N_field,Sizes_ghosted)
left_parity=get_left_parity(fields)
right_parity=get_right_parity(fields)
periodicity=get_periodicity(fields)
cache=FieldsCache(T,fields,N_copies)

DiscreteFileds{T,1+space_dimension,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,S,N_field*N_field}(
fields,discretization,idx,left_parity,right_parity,periodicity,cache
)

end 

function convert_field(type::Type{T},x::DiscreteFileds{T2,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,S,N_field2}) where {T,T2,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,S,N_field2}
    DiscreteFileds{T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,S,N_field2}(x.fields,x.discretization,x.index_structure,x.left_parity,x.right_parity,x.periodicity,convert_cache(T,x.cache))
end

function get_array(disc::DiscreteFileds{T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,S,N_field2}) where {T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,S,N_field2}
    zeros(T,(N_field,Sizes_ghosted...))
end


"""
    set_array!(array,fun,express::S,disc::DiscreteFileds) 

Set an array in place at the given field component specified by the espression with the function fun.
"""
function set_array!(array,fun,express::S,disc::DiscreteFileds{T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,S,N_field2}) where {S,T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,N_field2}
    #array=zeros(T,(N_field,Sizes_ghosted...))
    i=get_index(express,disc.fields)
    for I in disc.index_structure.interior
        #point=disc.discretization[I]
        point=disc.discretization.grid[I]
        array[i,I]= fun(point...)
    end

    boundary_condition!(array,disc)

    return array
end

"""
    set_array(fun,express::S,disc::DiscreteFileds) 

Set an array out of place at the given field component specified by the espression with the function fun.
"""
function set_array(fun,express::S,disc::DiscreteFileds{T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,S,N_field2}) where {S,T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,N_field2}
    array=get_array(disc)
    i=get_index(express,disc.fields)
    for I in disc.index_structure.interior
        #point=disc.discretization[I]
        point=disc.discretization.grid[I]
        array[i,I]= fun(point...)
    end
    boundary_condition!(array,disc)
    return array
end


function Base.getindex(phi::AbstractArray,index::Int,I::CartesianIndex{space_dimension},disc::DiscreteFileds{T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,S,N_field2}) where {space_dimension,S,T,total_dimensions,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,N_field2}
    
    for ndim in SOneTo(space_dimension)
        if I[ndim]==1
            e_i=disc.index_structure.unit_inidices_space[ndim]
            I_left_boundary= I+e_i
            I_end=end_index(ndim,I,Sizes_ghosted[ndim]) - e_i
            return disc.left_parity[ndim][index]*phi[index,I_left_boundary]+ disc.periodicity[ndim][index]*phi[index,I_end]

        elseif I[ndim]==Sizes_ghosted[ndim]
            e_i=disc.index_structure.unit_inidices_space[ndim]
            I_right_boundary=I-e_i
            I_first= first_index(ndim,I) +e_i
            return disc.right_parity[ndim][index]*phi[index,I_right_boundary] + disc.periodicity[ndim][index]*phi[index,I_first]
        else 
            return phi[index,I]
        end
    end

end 



 
function Base.getindex(phi::AbstractArray,I::CartesianIndex{space_dimension},disc::DiscreteFileds{T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,S,N_field2}) where {space_dimension,S,T,total_dimensions,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,N_field2}
    
    for ndim in SOneTo(space_dimension)
        if I[ndim]==1
            e_i=disc.index_structure.unit_inidices_space[ndim]
            I_left_boundary= I+e_i
            I_end=end_index(ndim,I,Sizes_ghosted[ndim]) - e_i
            
            return disc.left_parity[ndim] .* view(phi,:,I_left_boundary) .+ disc.periodicity[ndim] .* view(phi,:,I_end)

        elseif I[ndim]==Sizes_ghosted[ndim]
            e_i=disc.index_structure.unit_inidices_space[ndim]
            I_right_boundary=I-e_i
            I_first= first_index(ndim,I) +e_i
            
            return disc.right_parity[ndim] .* view(phi,:,I_right_boundary) + disc.periodicity[ndim] .* view(phi,:,I_first)
        else 
            return view(phi,:,I)
        end
    end

end


function boundary_condition!(phi,
    indexes::IndexStructure{dimension_tot,dimension_space,N_filed,TupleDim},
    left_parity,right_parity,periodicity) where {dimension_tot,dimension_space,N_filed,TupleDim}
    @inbounds for ndim in SOneTo(dimension_space) 
        ei=indexes.unit_inidices[ndim]

        
        @inbounds @simd for I in indexes.left_ghost[ndim]
            I_left_boundary= I+ei
            I_end=end_index(ndim+1,I,TupleDim[ndim])-ei
            phi[I]= left_parity[ndim][I[1]]*phi[I_left_boundary]+ periodicity[ndim][I[1]]*phi[I_end]
        end 
    
        @inbounds @simd for I in indexes.right_ghost[ndim]
            I_right_boundary=I- ei
            I_first=first_index(ndim+1,I)+ei
            phi[I] = right_parity[ndim][I[1]]*phi[I_right_boundary] + periodicity[ndim][I[1]]*phi[I_first]
            
        end 
    end 
end


boundary_condition!(phi,discrete_fields::DiscreteFileds{T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,S,N_field2}) where {T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,S,N_field2}= 
boundary_condition!(phi,discrete_fields.index_structure,discrete_fields.left_parity,discrete_fields.right_parity,discrete_fields.periodicity)


# this function is obsolate 
#function get_ghosted(discrete_fields::DiscreteFileds{T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,DXS,M,S,N_field2}) where {T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,DXS,M,S,N_field2}
#    grid=zeros(SVector{space_dimension,M},Sizes_ghosted)
#    non_ghosted=discrete_fields.discretization
#    
#    function fill_the_point(ndim,I)
#        if I[ndim]==1 
#            return entry=non_ghosted[ndim,I[ndim]]-DXS[ndim]
#        
#        elseif I[ndim]==Sizes_ghosted[ndim]
#        
#            return entry=non_ghosted[ndim,I[ndim]-2]+DXS[ndim]
#        else 
#
#            return entry=non_ghosted[ndim,I[ndim]-1]
#    
#        end 
#    end 
#   
#    for I in CartesianIndices(grid)
#        
#        point=ntuple(i->fill_the_point(i,I),Val{space_dimension}())
#        grid[I]= SVector{space_dimension,M}(point)
#    
#    end
#    return grid 
#end


@inbounds function basicupwinding(dphi,phi,t,discrete_fields::DiscreteFileds{T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Length,DXS,M,S,N_field2},matrix_equation!,params) where {T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Length,DXS,M,S,N_field2}
    #this is the discretization
    #@show t
    X=discrete_fields.discretization
    #then we have the index structure 
    idx=discrete_fields.index_structure
    # the interior 
    Interior=idx.interior
    # this are the unit vector (1,0,0,...)
    e_i=idx.unit_inidices_space
    #parity stuff 
    left_parity=discrete_fields.left_parity
    right_parity=discrete_fields.right_parity
    periodicity=discrete_fields.periodicity
    #chace matrices where to write the intermidiate results
    cache=discrete_fields.cache
    #this is a ntuple of matrices Ndim
    @inbounds A_i= cache.matrix[1]
    #this is just a vector that is need it for the source term 
    @inbounds Source= cache.vector[1][1]
    @inbounds nabla= cache.vector[2]
    @inbounds delta= cache.vector[3]
    @inbounds temp= cache.vector[4][1]

    #now let do some computation regaring the grid 
    deltax=Δx(X)
    inv_deltax=1 ./deltax
    #inv_deltax_2=inv_deltax *inv_deltax
    onehalf= one(T)/2
    two=2*one(T)

    #Loop over all the interior point I is a space iteration (x,y,z,....) the first index is not iterate
    @inbounds @simd for I in Interior
        # now is a N_field vector here i do not allocate nothing  
        ϕ=@views phi[:,I]
        dϕ=@views dphi[:,I]
        #here we compute all the matrices (matrices,Source,phi,t,x,params) 
        matrix_equation!(A_i,Source,ϕ,t,X[I],params)
        #now that we have the matrices in A_i as a tuple we loop over the dimension 
        @. dϕ =  - Source
        
        @inbounds @simd for dim in SOneTo{space_dimension}()
            #for each dimension i shift the index in the plus and minus dimension
            versor= e_i[dim]
            plus_index=I +versor
            minus_index=I -versor
            ϕ_plus = @views phi[:,plus_index]
            ϕ_minus= @views phi[:,minus_index]
            ∇_iϕ=nabla[dim]
            Δ_iϕ=delta[dim]
            # here we compute the derivative 
            invdx=inv_deltax[dim]
            @. ∇_iϕ = (ϕ_plus - ϕ_minus) * invdx * onehalf
            @. Δ_iϕ = (ϕ_plus + ϕ_minus -two *ϕ) * invdx
            #for n_field in SOneTo(N_field) 
            #∇_iϕ[n_field] = (ϕ_plus[n_field] - ϕ_minus[n_field]) * invdx * onehalf
            #Δ_iϕ[n_field] = (ϕ_plus[n_field] + ϕ_minus[n_field] -two *ϕ[n_field]) * invdx
            #end 
            #here we compute the flux
            upwindflux!(Δ_iϕ,∇_iϕ,A_i[dim],temp)
            #upwindflux!(Δ_iϕ,∇_iϕ,A_i[dim],temp)
            @inbounds @simd for n_field in SOneTo{N_field}()
                dϕ[n_field] -=  Δ_iϕ[n_field]
            end 
        end  

    end

    #before returning we  we apply the boundary_condition to the array are already decided 
    boundary_condition!(dphi,idx,left_parity, right_parity, periodicity)

    return nothing
end




@inline function upwindflux!(diffusion,nabla,derivativeMatrix,temp ) 
    jgemvavx!(temp,derivativeMatrix ,diffusion)

    jgemvavx!(diffusion, derivativeMatrix,temp,1,1)

    jgemvavx!(diffusion,derivativeMatrix,nabla,1,-1/4)
end

function cheb_coeff(j;N=3,f=abs)
    @assert isless(j,N)
    arg1(k) = π*(k+0.5)/N
    arg2(k,j) = π*j*(k+0.5)/N
    sum = 0
    [sum+=f(cos(arg1(k)))*cos(arg2(k,j)) for k in 0:N-1]
    return 2/N*sum
end


function cheb_pol(n,x)
    @assert isless(-1,n)
    if n == 0
        return one(x)
    elseif n==1
        return x
    else
        return 2*x*cheb_pol(n-1,x)-cheb_pol(n-2,x)
    end
end
#2*T2*T(2k-2)-T(2k-4)

function cheb_recursion!(T2k,T2,T2km2,T2km4)
    mul!(T2km4,T2,T2km2,2,-1)
    T2k .=  T2km4
end

@inline function abscoefficients(k)
    (-1)^(k+1)/(2k-1)/(2k+1)*2
end

@inline """
    mulladd_identiy!(dest::AbstractMatrix,in::AbstractMatrix,a,b)

Compute dest=in*a +b*I inplace 
"""
function mulladd_identiy!(dest::AbstractMatrix,in::AbstractMatrix,a,b)
    
    @inbounds for i in CartesianIndices(in)
        result=in[i]*a 
        if i[1]==i[2]
            result=result+b
        end
        dest[i]=result 
    end 
end

@inline """
    mulladd_identiy!(dest::AbstractMatrix,in::AbstractMatrix,a,b)

Compute in=in*a +b*I inplace 
"""
function mulladd_identiy!(in::AbstractMatrix,a,b)
    
    @inbounds for i in CartesianIndices(in)
        result=in[i]*a 
        if i[1]==i[2]
            result=result+b
        end
        in[i]=result 
    end 
end


function cheb_flux!(diffusion,A,N,T2=similar(A),T2k=similar(diffusion),T2km2=similar(diffusion),T2km4=similar(diffusion))
   
    mul!(T2,A,A)
    mulladd_identiy!(T2,2*one(eltype(A)),-one(eltype(A)))
    #T2 .= T2.*2 .- one(T2)
    
    mul!(T2km2,T2,diffusion)
    
    T2km4 .= diffusion
    
    fill!(T2k,zero(eltype(diffusion)))

    @. diffusion = T2km4 + abscoefficients(1)*T2km2
    count = 2

    while count<N
        cheb_recursion!(T2k,T2,T2km2,T2km4)
        @. diffusion = diffusion + abscoefficients(count)*T2km2
        count+=1
        T2km4 .= T2km2
        T2km2 .= T2k
    end
    diffusion .= 2/pi .*diffusion

    return nothing
end

function cheb_flux!(diffusion,nabla,A,N,T2=similar(A),T2k=similar(diffusion),T2km2=similar(diffusion),T2km4=similar(diffusion))
    cheb_flux!(diffusion,A,N,T2,T2k,T2km2,T2km4)
    mul!(diffusion,derivativeMatrix,nabla,1,-1/2)
end

@inbounds function chebupwinding(dphi,phi,t,discrete_fields::DiscreteFileds{T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Length,DXS,M,S,N_field2},matrix_equation!,params) where {T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Length,DXS,M,S,N_field2}
    #this is the discretization
    #@show t
    X=discrete_fields.discretization
    #then we have the index structure 
    idx=discrete_fields.index_structure
    # the interior 
    Interior=idx.interior
    # this are the unit vector (1,0,0,...)
    e_i=idx.unit_inidices_space
    #parity stuff 
    left_parity=discrete_fields.left_parity
    right_parity=discrete_fields.right_parity
    periodicity=discrete_fields.periodicity
    #chace matrices where to write the intermidiate results
    cache=discrete_fields.cache
    #this is a ntuple of matrices Ndim
    @inbounds A_i= cache.matrix[1]
    #this is just a vector that is need it for the source term 
    @inbounds Source= cache.vector[1][1]
    @inbounds nabla= cache.vector[2]
    @inbounds delta= cache.vector[3]

    #cache for intermidiate calculation 
    #@inbounds temp= cache.vector[4][1]
    @inbounds T2_mat=cache.matrix[2][1]
    @inbounds T2k_vec=cache.vector[4][1]
    @inbounds T2km2_vec=cache.vector[5][1]
    @inbounds T2km4_vec=cache.vector[6][1]

    #now let do some computation regaring the grid 
    deltax=Δx(X)
    inv_deltax=1 ./deltax
    #inv_deltax_2=inv_deltax *inv_deltax
    onehalf= one(T)/2
    two=2*one(T)

    #Loop over all the interior point I is a space iteration (x,y,z,....) the first index is not iterate
    @inbounds @simd for I in Interior
        # now is a N_field vector here i do not allocate nothing  
        ϕ=@views phi[:,I]
        dϕ=@views dphi[:,I]
        #here we compute all the matrices (matrices,Source,phi,t,x,params) 
        matrix_equation!(A_i,Source,ϕ,t,X[I],params)
        #now that we have the matrices in A_i as a tuple we loop over the dimension 
        @. dϕ =  - Source
        
        @inbounds @simd for dim in SOneTo{space_dimension}()
            #for each dimension i shift the index in the plus and minus dimension
            versor= e_i[dim]
            plus_index=I +versor
            minus_index=I -versor
            ϕ_plus = @views phi[:,plus_index]
            ϕ_minus= @views phi[:,minus_index]
            ∇_iϕ=nabla[dim]
            Δ_iϕ=delta[dim]
            # here we compute the derivative 
            invdx=inv_deltax[dim]
            @. ∇_iϕ = (ϕ_plus - ϕ_minus) * invdx * onehalf
            @. Δ_iϕ = (ϕ_plus + ϕ_minus -two *ϕ) * invdx
            #for n_field in SOneTo(N_field) 
            #∇_iϕ[n_field] = (ϕ_plus[n_field] - ϕ_minus[n_field]) * invdx * onehalf
            #Δ_iϕ[n_field] = (ϕ_plus[n_field] + ϕ_minus[n_field] -two *ϕ[n_field]) * invdx
            #end 
            #here we compute the flux
            #upwindflux!(Δ_iϕ,∇_iϕ,A_i[dim],temp)
            cheb_flux!(Δ_iϕ,∇_iϕ,A_i[dim],3,T2_mat,T2k_vec,T2km2_vec,T2km4_vec)
            #upwindflux!(Δ_iϕ,∇_iϕ,A_i[dim],temp)
            @inbounds @simd for n_field in SOneTo{N_field}()
                dϕ[n_field] -=  Δ_iϕ[n_field]
            end 
        end  

    end

    #before returning we  we apply the boundary_condition to the array are already decided 
    boundary_condition!(dphi,idx,left_parity, right_parity, periodicity)

    return nothing
end






function cheb_approx(x;N=3, f= abs)
    sum = zero(x)
    [sum+=cheb_coeff(j,N=N,f=f)*cheb_pol(j,x) for j in 0:N-1]
    return sum -one(x)*0.5*cheb_coeff(0,N=N,f=f)
end


@inline function upwindflux_Chebyschev!(diffusion,nabla,derivativeMatrix,temp ) 
    @turbo for i ∈ eachindex(diffusion)
        diffusioni = diffusion[i]
        for j ∈ eachindex(diffusion)
            diffusioni = -0.5*cheb_approx(derivativeMatrix[i,j]) * diffusion[j] + derivativeMatrix[i,j]*nabla[j]
        end
        diffusion[i] = diffusioni
    end
    
end


@inline function upwindfluxcomplex!(diffusion,nabla,derivativeMatrix,temp ) 
    
    mul!(temp,derivativeMatrix ,diffusion)

    mul!(diffusion, derivativeMatrix,temp,1,1)

    mul!(diffusion,derivativeMatrix,nabla,1,-1/4)
end




@inline function jgemvavx!(𝐲, 𝐀, 𝐱)
    @turbo for i ∈ eachindex(𝐲)
        𝐲i = zero(eltype(𝐲))
        for j ∈ eachindex(𝐱)
            𝐲i += 𝐀[i,j] * 𝐱[j]
        end
        𝐲[i] = 𝐲i
    end
end

@inline function jgemvavx!(𝐲, 𝐀, 𝐱,alpha,beta)
    @turbo for i ∈ eachindex(𝐲)
        𝐲i = beta*𝐲[i]
        for j ∈ eachindex(𝐱)
            𝐲i += alpha*𝐀[i,j] * 𝐱[j]
        end
        𝐲[i] = 𝐲i
    end
end


function problem(two_ideal_hydro_discrete::DiscreteFileds{T, total_dimensions, space_dimension, N_field, Sizes_ghosted, Sizes, Lengths, DXS, M, S, N_field2},matrix_eq,cs,init,tspan) where {T, total_dimensions, space_dimension, N_field, Sizes_ghosted, Sizes,Lengths, DXS, M, S, N_field2}
    #NewT=promote_type(eltype(init), typeof(cs))
    NewT=eltype(init)
    
    two_ideal_hydro_discrete=convert_field(NewT,two_ideal_hydro_discrete)
        ODEProblem((du,u,p,t)->basicupwinding(du,u,t,two_ideal_hydro_discrete,matrix_eq,cs),NewT.(init),tspan)
end

function oneshoot(two_ideal_hydro_discrete,ideal_matrix_equation_2d!,cs,phi,tspan,args...;kwargs...)
    prob=problem(two_ideal_hydro_discrete,ideal_matrix_equation_2d!,cs,phi,tspan)
    solve(prob,args...;kwargs...)
end



function isosurface(disc::DiscreteFileds{T, total_dimensions, space_dimension, N_field, Sizes_ghosted, Sizes,Lengths, DXS, M, S, N_field2},
    ideal_matrix_equation_2d!,cs,phi,tspan,express::Symbol,surf) where {T, total_dimensions, space_dimension, N_field, Sizes_ghosted, Sizes,Lengths, DXS, M, S, N_field2}
    
    prob=problem(disc,ideal_matrix_equation_2d!,cs,phi,tspan)
    
    integrator = init(prob,Tsit5();save_everystep=false)
    e_i=disc.index_structure.unit_inidices_space
    X=disc.discretization
    index=get_index(express,disc.fields)
#preallocation of a surface vector 
    surface_list=zeros(surface_crossing_point{T,typeof(tspan[1]),space_dimension+1,N_field},100000)
    max_size=length(surface_list)
    count=1
    max_of_index=zero(T)
    @inbounds for (uprev,tprev,u,t) in intervals(integrator)
        #@show t ,max_of_index  
        max_of_index=zero(T)
        @inbounds for I in disc.index_structure.interior
            ϕ=@views u[:,I]
            ϕprev=@views uprev[:,I]
            max_of_index=max(ϕprev[index],max_of_index)
            # time  
            #if ((ϕ[index]>surf)&&(ϕprev[index]<surf))||((ϕprev[index]>surf)&&(ϕ[index]<surf))
            if (ϕprev[index]>surf)&&(ϕ[index]<surf)
                if count<=max_size
                surface_list[count]=surface_crossing_point(X[tprev,I,Val{:SVector}()],X[t,I,Val{:SVector}()],SVector{N_field,T}(ntuple(i->ϕprev[i],Val{N_field}())),SVector{N_field,T}(ntuple(i->ϕ[i],Val{N_field}())))
                count=count +1 
                else 
                    push!(surface_list,surface_crossing_point(X[tprev,I,Val{:SVector}()],X[t,I,Val{:SVector}()],SVector{N_field,T}(ntuple(i->ϕprev[i],Val{N_field}())),SVector{N_field,T}(ntuple(i->ϕ[i],Val{N_field}()))))
                    count=count +1 
                end 
            end 
            #space (we check uprev)
            @inbounds for versor in e_i
                I_check_plus= I + versor
                I_check_minus= I - versor
                ϕ_plus = @views uprev[:,I_check_plus]
                ϕ_minus = @views uprev[:,I_check_minus]
                #if ((ϕ_plus[index]>surf)&&(ϕprev[index]<surf))||((ϕprev[index]>surf)&&(ϕ_plus[index]<surf))
                if (ϕprev[index]>surf)&&(ϕ_plus[index]<surf)
                    if count<=max_size
                        surface_list[count]=surface_crossing_point(
                            X[tprev,I,Val{:SVector}()],
                            X[tprev,I_check_plus,Val{:SVector}()],
                            SVector{N_field,T}(ntuple(i->ϕprev[i],Val{N_field}())),
                            SVector{N_field,T}(ntuple(i->ϕ_plus[i],Val{N_field}()))
                            )
                        count=count +1 
                    else 
                        push!(surface_list,surface_crossing_point(
                            X[tprev,I,Val{:SVector}()],
                            X[tprev,I_check_plus,Val{:SVector}()],
                            SVector{N_field,T}(ntuple(i->ϕprev[i],Val{N_field}())),
                            SVector{N_field,T}(ntuple(i->ϕ_plus[i],Val{N_field}())))
                            )
                        count=count +1
                    end 
                end 
                #if ((ϕ_minus[index]>surf)&&(ϕprev[index]<surf))||((ϕprev[index]>surf)&&(ϕ_minus[index]<surf))
                if ((ϕprev[index]>surf)&&(ϕ_minus[index]<surf))
                    if count<=max_size
                        surface_list[count]=surface_crossing_point(
                            X[tprev,I,Val{:SVector}()],
                            X[tprev,I_check_minus,Val{:SVector}()],
                            SVector{N_field,T}(ntuple(i->ϕprev[i],Val{N_field}())),
                            SVector{N_field,T}(ntuple(i->ϕ_minus[i],Val{N_field}()))
                            )
                        count=count +1  
                    else 
                        push!(surface_list,surface_crossing_point(
                            X[tprev,I,Val{:SVector}()],
                            X[tprev,I_check_minus,Val{:SVector}()],
                            SVector{N_field,T}(ntuple(i->ϕprev[i],Val{N_field}())),
                            SVector{N_field,T}(ntuple(i->ϕ_minus[i],Val{N_field}())))
                            )
                        count=count +1
                    end 
                end 
            end 


        end 
        
        if max_of_index<surf
            terminate!(integrator)
        end 
    end

    if length(surface_list)>count
        surface_list=surface_list[1:count-1]
    end  
    surface=linar_interpol.(surface_list,Ref(surf),Ref(express),Ref(disc))
   return (;surface, surface_list, surf ,integrator)
end





struct surface_crossing{S,T,N_dim,N_field}
    t_1::T
    t_2::T
    I_1::SVector{N_field,T}
    I_2::SVector{N_field,T}
    phi_1::SVector{N_field,S}
    phi_2::SVector{N_field,S}
end

function Base.zero(::Type{surface_crossing{S,T,N_dim,N_field}}) where {S,T,N_dim,N_field}
    surface_crossing(
        zero(T)
    ,zero(T),
    zero(CartesianIndex{N_dim}),
    zero(CartesianIndex{N_dim}),
    zeros(SVector{N_field,S}),
    zeros(SVector{N_field,S})  )
end

struct surface_crossing_point{S,T,N_dim,N_field}
    X_1::SVector{N_dim,T}
    X_2::SVector{N_dim,T}
    phi_1::SVector{N_field,S}
    phi_2::SVector{N_field,S}
end 


function linar_interpol(S::surface_crossing_point{M,T,N_dim,N_field},temp,express::Symbol,disc) where {M,T,N_dim,N_field}
    index=get_index(express,disc.fields)
    x1=S.X_1
    x2=S.X_2
    phi1=S.phi_1
    phi2=S.phi_2
    temp1=phi1[index]
    temp2=phi2[index]
    DXDT=(x1-x2)/(temp1-temp2)
    X_fo=x2+(temp-temp2)*DXDT
    xfox1=first(filter(x->!iszero(x),X_fo-x1))
    x2x1=first(filter(x->!iszero(x),x2-x1))

    #phi_fo=phi1 +(norm(X_fo-x1)/norm(x2-x1))*(phi2-phi1)
    phi_fo=phi1 +(xfox1/x2x1)*(phi2-phi1)

    surface_point(X_fo,phi_fo)
end

function Base.zero(::Type{surface_crossing_point{S,T,N_dim,N_field}}) where {S,T,N_dim,N_field}
surface_crossing_point(
    zeros(SVector{N_dim,T}),
    zeros(SVector{N_dim,T}),
zeros(SVector{N_field,S}),
zeros(SVector{N_field,S})  )
end

struct surface_point{S,T,N_dim,N_field}
    X::SVector{N_dim,T}
    phi::SVector{N_field,S}
end 

struct Surface_coodrinates{S,T,N_parm,N_dim,N_field}
    coordinates::SVector{N_parm,T}
    X::SVector{N_dim,T}
    phi::SVector{N_field,S}
end 


function Surface_coodrinates(point::surface_point{S,T,N_dim,N_field} ) where {S,T,N_dim,N_field}
    #here i pick up the first two entry of x 
    x_red=point.X[1:2]
    r=hypot(x_red...)
    theta=atan(x_red...)
    if N_dim <3
        Surface_coodrinates(SVector{N_dim-1}(theta),point.X,point.phi
        )
    else 
        Surface_coodrinates(SVector{N_dim-1}((theta,point.X[3:end]...)),
        point.X,point.phi
        )
    end 
end 


function Surface_coodrinates(point::surface_point{S,T,N_dim,N_field},surface_condition ) where {S,T,N_dim,N_field}
    #here i pick up the first two entry of x 
    #x_red=point.X[1:2]
    #r=hypot(x_red...)
    #theta=atan(x_red...)
    
    coordianates=surface_condition(point.X...)
    
    Surface_coodrinates(coordianates,point.X,point.phi
        )
end


struct Surface{S,T,N_dim,N_field}
    points::Vector{surface_point{S,T,N_dim,N_field}}
end 

struct Chart{S,T,N_parm,N_dim,N_field}
    points::Vector{Surface_coodrinates{S,T,N_parm,N_dim,N_field}}
end 


function Chart(Sur::Surface{S,T,N_dim,N_field}) where {S,T,N_dim,N_field}
    Chart{S,T,N_dim-1,N_dim,N_field}(Surface_coodrinates.(Sur.points))
end

function Chart(Sur::Surface{S,T,N_dim,N_field},surface_condition) where {S,T,N_dim,N_field}
    Chart{S,T,N_dim-1,N_dim,N_field}(Surface_coodrinates.(Sur.points,surface_condition))
end


struct DifferentibleChart{S,T,N_parm,N_dim,N_field,leftbound,rightbound,A,B}
    X::A
    phi::B
end

function DifferentibleChart(Sur::Chart{S,T,N_parm,N_dim,N_field},decimation::Int,cheb_tuple::NTuple{N_parm,Int}) where  {S,T,N_parm,N_dim,N_field}
    
    resample=StructArray(StatsBase.sample(Sur.points, length(Sur.points)÷decimation; replace=false))

    X=resample.X
    alpha=resample.coordinates

    phi=resample.phi
    left=ntuple(i->extrema(getindex.(resample.coordinates,(i)))[1],Val{N_parm}())

    right=ntuple(i->extrema(getindex.(resample.coordinates,(i)))[2],Val{N_parm}())
    
    X_interp=ntuple(i->RadialBasis(alpha, getindex.(X,(i)),left,right,rad=thinplateRadial()),Val{N_dim}())

    phi_interp=ntuple(i->RadialBasis(alpha, getindex.(phi,(i)),left,right,rad=thinplateRadial()),Val{N_dim}())

    function X_threaded(alpha)
        SVector{N_dim,T}(ntuple(i->X_interp[i](alpha) ,Val{N_dim}()))
    end 

    function phi_threaded(alpha)
        SVector{N_dim,S}(ntuple(i->phi_interp[i](alpha) ,Val{N_dim}()))
    end 

    re_alpha = chebpoints(cheb_tuple, left, right)
    #re_X = chebinterp(X_threaded.(re_alpha), left, right)

    re_X = chebregression(re_alpha,X_threaded.(re_alpha...), left, right,70)

    #re_phi = chebinterp(phi_threaded.(re_alpha), left, right)

    re_phi = chebregression(re_alpha,phi_threaded.(re_alpha...), left, right,70)

    DifferentibleChart{S,T,N_parm,N_dim,N_field,left,right,typeof(re_X),typeof(re_phi)}(re_X,re_phi)
end 

function left_bound(dchart::DifferentibleChart{S,T,N_parm,N_dim,N_field,leftbound,rightbound,A,B}) where {S,T,N_parm,N_dim,N_field,leftbound,rightbound,A,B}
    return leftbound
end

function right_bound(dchart::DifferentibleChart{S,T,N_parm,N_dim,N_field,leftbound,rightbound,A,B}) where {S,T,N_parm,N_dim,N_field,leftbound,rightbound,A,B}
    return rightbound
end

struct FreezeOutResult{A,B}
    x::A
    fields::B
end 

Base.iterate(S::FreezeOutResult) = (S.x, Val(:x))
Base.iterate(S::FreezeOutResult, ::Val{:x}) = (S.fields, Val(:done))
Base.iterate(S::FreezeOutResult, ::Val{:done}) = nothing


function _radial_basisinterpolate(surf::AbstractVector{Surface_coodrinates{S,T,N_parm,N_dim,N_field}};sort_index=3) where {S,T,N_parm,N_dim,N_field}
    
    sortedcha=StructArray(sort(surf,by=v->getindex(v.coordinates,sort_index)))

    totalpoint=length(sortedcha )

    total_left=ntuple(i->extrema(getindex.(sortedcha.coordinates,(i)))[1],Val{N_parm}())

    total_right=ntuple(i->extrema(getindex.(sortedcha.coordinates,(i)))[2],Val{N_parm}())

    X_interp=ntuple(i->RadialBasis(sortedcha.coordinates, getindex.(sortedcha.X,(i)),total_left,total_right,rad=thinplateRadial()),Val{N_dim}())

    phi_interp=ntuple(i->RadialBasis(sortedcha.coordinates, getindex.(sortedcha.phi,(i)),total_left,total_right,rad=thinplateRadial()),Val{N_field}())

    #(;X_interp , phi_interp,  total_left ,total_right,totalpoint)

    function X(x)
        SVector{N_dim}(ntuple(i->X_interp[i](x),Val{N_dim}()))
    end

    function phi(x)
        SVector{N_field}(ntuple(i->phi_interp[i](x),Val{N_field}()))
    end

    (coord=BoundedFunction(X,total_left,total_right),fields=BoundedFunction(phi,total_left,total_right))
end


function radial_basisinterpolate(surf::Chart{S,T,N_parm,N_dim,N_field};baches=2001,sort_index=3) where {S,T,N_parm,N_dim,N_field}

        #here i sort the vector on the dimension select by sort_index 
        sortedcha=sort(surf.points,by=v->getindex(v.coordinates,sort_index))
        
        #first we check the size of the vector 
        len=length(sortedcha)

        #then we compute the number of baches 
        if len >baches
            # here we have at least 2 baches 
            nbaches=div(len,baches) +1 
                
            x,phi=_radial_basisinterpolate(view(sortedcha,1:baches),sort_index=sort_index)
            
            X=PieceWiseFunction(x)
            Phi=PieceWiseFunction(phi)
            for i in 2:nbaches
                #@show i
                if i==nbaches
                    final_lenght=len
                else 
                    final_lenght = (baches*(i))
                end 
                #@show final_lenght
                x,phi=_radial_basisinterpolate(view(sortedcha,(baches*(i-1)+1):final_lenght),sort_index=sort_index)
                
                 
                X=PieceWiseFunction(X,x)
                Phi=PieceWiseFunction(Phi,phi)

            end 
        else 
            x,phi=_radial_basisinterpolate(sortedcha,sort_index=sort_index)

            X=PieceWiseFunction(x)
            Phi=PieceWiseFunction(phi)
        end 

    return FreezeOutResult(X,Phi) #(coord=X,fields=Phi)
end 


function spline_interpolation(X::FreezeOutResult{A,B};ndim_tuple=100) where {A<:PieceWiseFunction,B<:PieceWiseFunction}
    
    #(coord= SplineInterp(x[:coord],ndim_tuple),fields=SplineInterp(x[:fields],ndim_tuple))
    x,phi=X
    FreezeOutResult(SplineInterp(x,ndim_tuple),SplineInterp(phi,ndim_tuple))
end 


function freezeout_interpolation(surf::Chart{S,T,N_parm,N_dim,N_field};baches=5001,sort_index=2,ndim_tuple=100) where {S,T,N_parm,N_dim,N_field}
     
    spline_interpolation(radial_basisinterpolate(surf,baches=baches,sort_index=sort_index),ndim_tuple=ndim_tuple)
end 

