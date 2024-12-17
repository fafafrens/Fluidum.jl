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

```julia
using Fluidum
eos = Heavy_Quark()
eos_HQ = HadronResonaceGas()
params=Fluidum.FluidProperties(eos,QGPViscosity(0.2,0.2),SimpleBulkViscosity(0.1,15.0),HQdiffusion(0.2,1.5))
fo=Fluidum.runFluidum_fo(Fluidum.Step_Intial_Condition(31.57,4),Fluidum.pQCD_Initial_Condition(1,70.,0.463),Fluidum.Trento_Intial_Condition(1),params,eos_HQ,0.4)
obs = Fluidum.Observables(fo,0.200,params,0.156)
```
