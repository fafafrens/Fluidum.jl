abstract type TransportProperty end

abstract type ShearViscosity <:TransportProperty end 

abstract type BulkViscosity <:TransportProperty end 

abstract type Diffusion <:TransportProperty end 



struct FluidProperties{A<:EquationOfState,B<:ShearViscosity,C<:BulkViscosity,D<:Diffusion}
    eos::A
    shear::B
    bulk::C
    diffusion::D
end



thermodynamic(T,x::FluidProperties{A,B,C,D}) where {A,B,C,D}= thermodynamic(T,x.eos)
thermodynamic(T,α,x::FluidProperties{A,B,C,D}) where {A,B,C,D}= thermodynamic(T,α,x.eos)

pressure(T,wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= pressure(T,wal.eos)

pressure_derivative(T,::Val{1},wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= pressure_derivative(T,Val{1}(),wal.eos)
pressure_derivative(T,::Val{2},wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= pressure_derivative(T,Val{2}(),wal.eos)

pressure(T,α,wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= pressure(T,α,wal.eos)

pressure_derivative(T,α,::Val{1},::Val{0},wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= pressure_derivative(T,α,Val{1}(),Val{0}(),wal.eos)
pressure_derivative(T,α,::Val{0},::Val{1},wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= pressure_derivative(T,α,Val{0}(),Val{1}(),wal.eos)
pressure_derivative(T,α,::Val{2},::Val{0},wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= pressure_derivative(T,α,Val{2}(),Val{0}(),wal.eos)
pressure_derivative(T,α,::Val{1},::Val{1},wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= pressure_derivative(T,α,Val{1}(),Val{1}(),wal.eos)
pressure_derivative(T,α,::Val{0},::Val{2},wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= pressure_derivative(T,α,Val{0}(),Val{2}(),wal.eos)
 
energy_density(T,α,wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= energy_density(T,α,wal.eos)

viscosity(T,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=viscosity(T,thermodynamic(T,x.eos),x.shear)
viscosity(T,α,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=viscosity(T,α,thermodynamic(T,α,x.eos),x.shear)

τ_shear(T,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=τ_shear(T,thermodynamic(T,x.eos),x.shear)
τ_shear(T,α,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=τ_shear(T,α,thermodynamic(T,α,x.eos),x.shear)


viscosity(T,x::EquationOfState,y::ShearViscosity)=viscosity(T,thermodynamic(T,x),y)
viscosity(T,α,x::EquationOfState,y::ShearViscosity)=viscosity(T,α,thermodynamic(T,α,x),y)

τ_shear(T,x::EquationOfState,y::ShearViscosity)=τ_shear(T,thermodynamic(T,x),y)
τ_shear(T,α,x::EquationOfState,y::ShearViscosity)=τ_shear(T,α,thermodynamic(T,α,x),y)



bulk_viscosity(T,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=bulk_viscosity(T,thermodynamic(T,x.eos),x.bulk)
bulk_viscosity(T,α,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=bulk_viscosity(T,α,thermodynamic(T,α,x.eos),x.bulk)

τ_bulk(T,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=τ_bulk(T,thermodynamic(T,x.eos),x.bulk)
τ_bulk(T,α,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=τ_bulk(T,α,thermodynamic(T,α,x.eos),x.bulk)


bulk_viscosity(T,x::EquationOfState,y::BulkViscosity)=bulk_viscosity(T,thermodynamic(T,x),y)
bulk_viscosity(T,α,x::EquationOfState,y::BulkViscosity)=bulk_viscosity(T,α,thermodynamic(T,α,x),y)

τ_bulk(T,x::EquationOfState,y::BulkViscosity)=τ_bulk(T,thermodynamic(T,x),y)
τ_bulk(T,α,x::EquationOfState,y::BulkViscosity)=τ_bulk(T,α,thermodynamic(T,α,x),y)



diffusion(T,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=diffusion(T,thermodynamic(T,x.eos),x.diffusion)
diffusion(T,α,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=diffusion(T,α,thermodynamic(T,α,x.eos),x.diffusion)

τ_diffusion(T,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=τ_diffusion(T,thermodynamic(T,x.eos),x.diffusion)
τ_diffusion(T,α,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=τ_diffusion(T,α,thermodynamic(T,α,x.eos),x.diffusion)


diffusion(T,x::EquationOfState,y::Diffusion)=diffusion(T,thermodynamic(T,x.eos),y)
diffusion(T,α,x::EquationOfState,y::Diffusion)=diffusion(T,α,thermodynamic(T,α,x.eos),y)

τ_diffusion(T,x::EquationOfState,y::Diffusion)=τ_diffusion(T,thermodynamic(T,x.eos),y)
τ_diffusion(T,α,x::EquationOfState,y::Diffusion)=τ_diffusion(T,α,thermodynamic(T,α,x.eos),y)




