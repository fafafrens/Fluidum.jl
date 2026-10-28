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

## Example 1+1 dimension
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
    
disc_fields = DiscreteFields(oned_visc_hydro,disc,Float64) 
#define a example function for the temperature
function temperature(r)
       0.4*1/(exp(r/7)+1 )+0.0001
end
#set the temperature field
phi=set_array((x)->temperature(x),:temperature,disc_fields); 
#define a time range for the evolution
tspan = (0.2,1000)
#evolve until freeze-out temperature Tfo
freeze_out_routine(disc_fields,Fluidum.matrix1d_visc!,params,phi,tspan;Tfo=0.1565)
```

## Example 2+1 dimension 
Instead you can also solve in 2 dimension 
```julia
using Fluidum
#set up the fluid poroperties 
fluidpropery=FluidProperties(FluiduMEoS(),SimpleShearViscosity(0.1,.1),ZeroBulkViscosity(),ZeroDiffusion())

# the convention here are T, ux, uy, piyy, pizz, pixy, piB this has to match with the matrix defined
twod_visc_hydro=Fields(
NDField((:ghost,:ghost),(:ghost,:ghost),:temperature),
NDField((:ghost,:ghost),(:ghost,:ghost),:ux),
NDField((:ghost,:ghost),(:ghost,:ghost),:uy),
NDField((:ghost,:ghost),(:ghost,:ghost),:piyy),
NDField((:ghost,:ghost),(:ghost,:ghost),:pizz),
NDField((:ghost,:ghost),(:ghost,:ghost),:pixy),
NDField((:ghost,:ghost),(:ghost,:ghost),:piB)
)

#we define a 2 cartesian grid form -25 to 25 50 point each dimension 
discretization=CartesianDiscretization(Fluidum.SymmetricInterval(50,25.),Fluidum.SymmetricInterval(50,25.))

# we prepare the field with the discretization
twod_visc_hydro_discrete=DiscreteFields(twod_visc_hydro,discretization,Float64)

 #we define some random intial condition 
function temperature(r)
       0.4(1+0.3rand())/(exp(abs(r)-10)+1 )+0.01
end
#we set the array corresponding to the temperature 
phi=set_array((x,y)->temperature(hypot(x,y)),:temperature,twod_visc_hydro_discrete);

using Plots
#here to undertand whath we are doing 
plot(phi[1,:,25])

#this is the time span 
tspan=(0.4,20)


# this create the ODEProblem  and feed it diffenential equations 
res=oneshoot(twod_visc_hydro_discrete,Fluidum.matrix2d_visc!,fluidpropery,phi,tspan)

#plot the solution 
using Plots
plot(phi[1,:,25])
plot!(res[end][1,:,25])
```

## Threads
Start Julia with several threads and the solver updates the cells of the grid in parallel:
```sh
julia -t 8
```
Nothing in the code changes, and the result is bit for bit the same as with one thread.

## Maximum time step
`oneshoot` limits the internal time step to `dtmax = 0.01` fm/c. Pass a smaller value as a keyword,
for example in the 2+1 dimensional example above:
```julia
res = oneshoot(twod_visc_hydro_discrete, Fluidum.matrix2d_visc!, fluidpropery, phi, tspan; dtmax = 0.005)
```
or set it for every solve of a run with the environment variable `FLUIDUM_DTMAX=0.005`.

## Medium and charm together in 2+1 dimension
`two_way_problem` evolves the viscous medium and the charm quarks (density, diffusion current and stress)
as one system, here from a fluctuating initial state of random hot spots. The charm field is the fugacity α.
It can be set directly or from a charm density, and the charm density follows from T and α. The charm stress
starts with less longitudinal than transverse momentum. The plot shows the fields T and α and the charm density
at the start and at the end.
```julia
using Fluidum, OrdinaryDiffEq, Random, Plots

# medium: lattice EoS, shear and bulk viscosity; charm: D_s T = 0.116 for a quark of mass 1.5 GeV
params = FluidProperties(FluiduMEoS(), QGPViscosity(0.1, 0.2), SimpleBulkViscosity(0.083, 15.0),
                         ConstDiffusion(0.116, 1.5))
quarks = charm_quark_list()   # one charm quark

# a fluctuating initial state: 15 Gaussian hot spots at random positions (fm)
Random.seed!(2)
spots = [(2.5 .* randn(2), 0.5 + rand()) for _ in 1:15]
s(x, y) = sum(w * exp(-((x - c[1])^2 + (y - c[2])^2) / 1.2) for (c, w) in spots)
T0(x, y) = 0.06 + 0.16 * cbrt(s(x, y))     # GeV

# 14 fields: T, ux, uy, piyy, pizz, pixy, piB, α, nux, nuy, piyy_Q, pizz_Q, pixy_Q, piB_Q
L, N = 10.0, 60
d = DiscreteFields(two_way_2d(), CartesianDiscretization(Fluidum.SymmetricInterval(N, L),
                                                         Fluidum.SymmetricInterval(N, L)), Float64)
phi = set_array(T0, :temperature, d)

# the charm fugacity α, set directly: α = 0 is charm in chemical equilibrium with the medium
α0(x, y) = 0.0
# or from a charm density n0, here produced in the same hot spots; for one quark n is proportional to e^α
# n0(x, y) = 0.05 * s(x, y) + 1e-8   # charm per fm^3
# α0(x, y) = log(n0(x, y) / Fluidum.charm_eos(T0(x, y), 0.0, quarks).n)
set_array!(phi, α0, :α, d)

# the charm stress (second moment), in units of the charm pressure P0 = n T: charm produced with little
# longitudinal momentum, pizz_Q = -a P0 and piyy_Q = a P0 / 2 (pixx_Q follows from tracelessness);
# pixy_Q and piB_Q stay zero, as does the medium's stress
a = 0.5
P0(x, y) = Fluidum.charm_eos(T0(x, y), α0(x, y), quarks).n * T0(x, y)
set_array!(phi, (x, y) -> -a * P0(x, y), :pizz_Q, d)
set_array!(phi, (x, y) -> a / 2 * P0(x, y), :piyy_Q, d)

# medium and charm evolved together; the default ColdEdge decouples charm below T = 0.156 GeV
prob = two_way_problem(d, params, phi, (0.4, 6.0); species = quarks)
sol = solve(prob, Tsit5(); saveat = [0.4, 6.0])

# the fields T and α, and the charm density computed from them, at the start and at the end
xs = range(-L + L / N, L - L / N; length = N)   # the cell centres
T(u) = u[1, 2:end-1, 2:end-1]
α(u) = u[8, 2:end-1, 2:end-1]
n(u) = [Fluidum.charm_eos(u[1, i, j], u[8, i, j], quarks).n for i in 2:N+1, j in 2:N+1]
plot(heatmap(xs, xs, T(sol.u[1])', title = "T (GeV), τ = 0.4 fm/c"),
     heatmap(xs, xs, T(sol.u[end])', title = "T (GeV), τ = 6 fm/c"),
     heatmap(xs, xs, α(sol.u[1])', title = "α, τ = 0.4 fm/c"),
     heatmap(xs, xs, α(sol.u[end])', title = "α, τ = 6 fm/c"),
     heatmap(xs, xs, n(sol.u[1])', title = "charm density (fm⁻³), τ = 0.4 fm/c"),
     heatmap(xs, xs, n(sol.u[end])', title = "charm density (fm⁻³), τ = 6 fm/c"),
     layout = (3, 2), size = (900, 1150), aspect_ratio = 1, xlabel = "x (fm)", ylabel = "y (fm)",
     left_margin = 5Plots.mm)
```

## Citing 
When using Fluidum.jl for research, teaching or similar, please cite our work, see [CITATIONS.bib](CITATIONS.bib).
