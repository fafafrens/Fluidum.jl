

struct CartesianDiscretization{N_dim,Sizes,Lengths,DXS,T}
    grid::Array{NTuple{N_dim,T},N_dim}
end 


function CartesianDiscretization(Sizes::NTuple{N,T},Lengths::NTuple{N,NTuple{2,S}}) where {N,T,S}
    
    ranges= Base.OneTo.(2 .+Sizes)
    lefts=Tuple( Lengths[i][1] for i in SOneTo{N}())
    rights=Tuple( Lengths[i][2] for i in SOneTo{N}())
    Dxs=Tuple( (Lengths[i][2]- Lengths[i][1]) / Sizes[i] for i in SOneTo{N}())

    grid=[ Dxs.* (Tuple(I) .-1) .+ lefts .- Dxs ./2   for I in CartesianIndices(ranges)]
    newL=Tuple( Tuple( convert(eltype(eltype(grid)),j) for j in i ) for i in Lengths)
    CartesianDiscretization{N,Sizes,newL,Dxs,eltype(eltype(grid))}(grid)
end 


@inline """
    CartesianDiscretization(oned_size::Int,oned_lenght::NTuple{2,T}) where {T}

Construct a one-dimensional grid of size oned_size and oned_lengh as extrama with ghosted cell at the boundary
"""
function CartesianDiscretization(oned_size::Int,left::X,right::Y) where {X<:Number,Y<:Number}

    CartesianDiscretization((oned_size,),(promote(left,right),))
end
@inline CartesianDiscretization(x::CartesianDiscretization{N_dim1,Sizes1,Lengths1,DXS1,T1}) where {N_dim1,Sizes1,Lengths1,DXS1,T1} =x

@inline function CartesianDiscretization(::CartesianDiscretization{N_dim1,Sizes1,Lengths1,DXS1,T1},::CartesianDiscretization{N_dim2,Sizes2,Lengths2,DXS2,T2}) where {N_dim1,Sizes1,Lengths1,DXS1,T1,N_dim2,Sizes2,Lengths2,DXS2,T2}
   
    
    CartesianDiscretization((Sizes1...,Sizes2...),(Lengths1...,Lengths2...))
end



"""
    CartesianDiscretization(x...)

Compose multiple cartesian grid as a tensor product.
"""
function CartesianDiscretization(x...)
    @show x
    CartesianDiscretization(
        CartesianDiscretization((first(x))
        ),CartesianDiscretization(Base.tail(x)...))
end
@inline Base.getindex(S::CartesianDiscretization{N_dim,Sizes,Lengths,DXS,T},i::Int) where {N_dim,Sizes,Lengths,DXS,T} = S.grid[i]
#@inline Base.firstindex(S::Interval) = 1
#@inline Base.lastindex(S::Interval{N,DX,T}) where {N,DX,T}=N 
@inline """
    Base.getindex(S::CartesianDiscretization{N_dim,Sizes,Lengths,DXS,T},i::CartesianIndex{N_dim}) where {N_dim,Sizes,Lengths,DXS,T}

Return a NTuple of the point corrsponding of the index i
"""
Base.getindex(S::CartesianDiscretization{N_dim,Sizes,Lengths,DXS,T},i::CartesianIndex{N_dim}) where {N_dim,Sizes,Lengths,DXS,T} =  S.grid[i]
@inline """
    Base.getindex(S::CartesianDiscretization{N_dim,Sizes,Lengths,DXS,T},i::CartesianIndex{N_dim},::Val{:SVector}) where {N_dim,Sizes,Lengths,DXS,T}

Return a SVector of the point corrsponding of the index i
"""
Base.getindex(S::CartesianDiscretization{N_dim,Sizes,Lengths,DXS,T},i::CartesianIndex{N_dim},::Val{:SVector}) where {N_dim,Sizes,Lengths,DXS,T} = SVector{N_dim,T}(S.grid[i])

@inline """
    Base.getindex(S::CartesianDiscretization{N_dim,Sizes,Lengths,DXS,T},t::M,i::CartesianIndex{N_dim},::Val{:SVector}) where {N_dim,Sizes,Lengths,DXS,T,M}

Return a SVector of the point corrsponding of the index i adding in front the value t 
"""
Base.getindex(S::CartesianDiscretization{N_dim,Sizes,Lengths,DXS,T},t::M,i::CartesianIndex{N_dim},::Val{:SVector}) where {N_dim,Sizes,Lengths,DXS,T,M} = SVector{N_dim+1,promote_type(T,M)}(t,S.grid[i]...)


@inline Base.getindex(S::CartesianDiscretization{N_dim,Sizes,Lengths,DXS,T},i::NTuple{N_dim,Int}) where {N_dim,Sizes,Lengths,DXS,T} = S.grid[i...]
@inline Base.getindex(S::CartesianDiscretization{N_dim,Sizes,Lengths,DXS,T},i::NTuple{N_dim,Int},::Val{:SVector}) where {N_dim,Sizes,Lengths,DXS,T} = SVector{N_dim,T}(S.grid[i...])

@inline """
    Δx(::CartesianDiscretization{N_dim,Sizes,Lengths,DXS,T}) where {N_dim,Sizes,Lengths,DXS,T}

Return the Δ of the cell in each dimension 
"""
Δx(::CartesianDiscretization{N_dim,Sizes,Lengths,DXS,T}) where {N_dim,Sizes,Lengths,DXS,T}=DXS


"""
    SymmetricInterval(oned_size::Int,oned_lenght::X) where{X<:Number}

Construct 1-D grid of size oned_size and oned_lenght as size simmetrically respect to the origin.
"""
SymmetricInterval(oned_size::Int,oned_lenght::X) where{X<:Number} =CartesianDiscretization(oned_size,-oned_lenght,oned_lenght)
"""
    PeriodicInterval(oned_size::Int)

Construct a 1-D grid from (0 ,2 π) with oned_size cells
"""
PeriodicInterval(oned_size::Int)=CartesianDiscretization(oned_size,0,2*pi)
"""
    OriginInterval(oned_size::Int,l::X) where{X<:Number}

Construct a 1-D grid from (0 , l) with oned_size cells
"""
OriginInterval(oned_size::Int,l::X) where{X<:Number}=CartesianDiscretization(oned_size,zero(l),l)

