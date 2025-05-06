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
thermodynamic(T,mu,x::FluidProperties{A,B,C,D}) where {A,B,C,D}= thermodynamic(T,mu,x.eos)

pressure(T,wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= pressure(T,wal.eos)

pressure_derivative(T,::Val{1},wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= pressure_derivative(T,Val{1}(),wal.eos)
pressure_derivative(T,::Val{2},wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= pressure_derivative(T,Val{2}(),wal.eos)

pressure(T,mu,wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= pressure(T,mu,wal.eos)

pressure_derivative(T,mu,::Val{1},::Val{0},wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= pressure_derivative(T,mu,Val{1}(),Val{0}(),wal.eos)
pressure_derivative(T,mu,::Val{0},::Val{1},wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= pressure_derivative(T,mu,Val{0}(),Val{1}(),wal.eos)
pressure_derivative(T,mu,::Val{2},::Val{0},wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= pressure_derivative(T,mu,Val{2}(),Val{0}(),wal.eos)
pressure_derivative(T,mu,::Val{1},::Val{1},wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= pressure_derivative(T,mu,Val{1}(),Val{1}(),wal.eos)
pressure_derivative(T,mu,::Val{0},::Val{2},wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= pressure_derivative(T,mu,Val{0}(),Val{2}(),wal.eos)
 
energy_density(T,mu,wal::FluidProperties{A,B,C,D}) where {A,B,C,D}= energy_density(T,mu,wal.eos)

viscosity(T,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=viscosity(T,thermodynamic(T,x.eos),x.shear)
viscosity(T,mu,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=viscosity(T,mu,thermodynamic(T,mu,x.eos),x.shear)

τ_shear(T,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=τ_shear(T,thermodynamic(T,x.eos),x.shear)
τ_shear(T,mu,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=τ_shear(T,mu,thermodynamic(T,mu,x.eos),x.shear)


viscosity(T,x::EquationOfState,y::ShearViscosity)=viscosity(T,thermodynamic(T,x),y)
viscosity(T,mu,x::EquationOfState,y::ShearViscosity)=viscosity(T,mu,thermodynamic(T,mu,x),y)

τ_shear(T,x::EquationOfState,y::ShearViscosity)=τ_shear(T,thermodynamic(T,x),y)
τ_shear(T,mu,x::EquationOfState,y::ShearViscosity)=τ_shear(T,mu,thermodynamic(T,mu,x),y)



bulk_viscosity(T,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=bulk_viscosity(T,thermodynamic(T,x.eos),x.bulk)
bulk_viscosity(T,mu,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=bulk_viscosity(T,mu,thermodynamic(T,mu,x.eos),x.bulk)

τ_bulk(T,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=τ_bulk(T,thermodynamic(T,x.eos),x.bulk)
τ_bulk(T,mu,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=τ_bulk(T,mu,thermodynamic(T,mu,x.eos),x.bulk)


bulk_viscosity(T,x::EquationOfState,y::BulkViscosity)=bulk_viscosity(T,thermodynamic(T,x),y)
bulk_viscosity(T,mu,x::EquationOfState,y::BulkViscosity)=bulk_viscosity(T,mu,thermodynamic(T,mu,x),y)

τ_bulk(T,x::EquationOfState,y::BulkViscosity)=τ_bulk(T,thermodynamic(T,x),y)
τ_bulk(T,mu,x::EquationOfState,y::BulkViscosity)=τ_bulk(T,mu,thermodynamic(T,mu,x),y)



diffusion(T,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=diffusion(T,thermodynamic(T,x.eos),x.diffusion)
diffusion(T,mu,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=diffusion(T,mu,thermodynamic(T,mu,x.eos),x.diffusion)

τ_diffusion(T,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=τ_diffusion(T,thermodynamic(T,x.eos),x.diffusion)
τ_diffusion(T,mu,x::FluidProperties{A,B,C,D}) where {A,B,C,D}=τ_diffusion(T,mu,thermodynamic(T,mu,x.eos),x.diffusion)


diffusion(T,x::EquationOfState,y::Diffusion)=diffusion(T,thermodynamic(T,x.eos),y)
diffusion(T,mu,x::EquationOfState,y::Diffusion)=diffusion(T,mu,thermodynamic(T,mu,x.eos),y)

τ_diffusion(T,x::EquationOfState,y::Diffusion)=τ_diffusion(T,thermodynamic(T,x.eos),y)
τ_diffusion(T,mu,x::EquationOfState,y::Diffusion)=τ_diffusion(T,mu,thermodynamic(T,mu,x.eos),y)




