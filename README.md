# Fluidum

[![Build Status](https://github.com/fafafrens/Fluidum.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fafafrens/Fluidum.jl/actions/workflows/CI.yml?query=branch%3Amain)

To install the package and its dependencies one can use 
```julia
import Pkg
Pkg.add(url = "https://github.com/fafafrens/Fluidum.jl")
```
After this one can use the package simply by using 
```julia
using Fluidum
```
## Example
A simple example of running heavy quark simulations 
(actually is not simple at all since it is not understandable: what Observable does?
what runFluidum ? 
Why we need 2 equations of state ? 
)

```julia
using Fluidum
#define equation of state and transport coefficients
eos = Heavy_Quark()
viscosity=QGPViscosity(0.2,0.2)
bulk=SimpleBulkViscosity(0.1,15.0)
diffusion=ZeroDiffusion()
#collect eos and transport coefficient in the FluidProperties struct
params=Fluidum.FluidProperties(eos,viscosity,bulk,diffusion)
#define a radial grid
gridpoints=100
rmax=30
disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
#define the name and number of fields, with appropriate parity
oned_visc_hydro = oned_visc_hydro=Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:piphiphi),
    NDField((:even,),(:ghost,),:pietaeta),
    NDField((:even,),(:ghost,),:piB)
    )
    
disc_fields = DiscreteFileds(oned_visc_hydro,disc,Float64) 
#define a example function for the temperature
function temperature(r)
       0.4*1/(exp(r/7)+1 )+0.0001
end
#set the temperature field
phi=set_array((x)->temperature(x),:temperature,disc_fields); 
#define a time range for the evolution
tspan = (0.2,1000)
#evolve until freeze-out temperature Tfo
freeze_out_routine(disc_fields,Fluidum.matrxi1d_visc!,params,phi,tspan;Tfo=0.1565)
```
