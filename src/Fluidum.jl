module Fluidum

using ForwardDiff:Dual,value,partials
using Bessels
using ForwardDiff
using DelimitedFiles: readdlm
using DelimitedFiles
using StructArrays
using IntervalSets
using MuladdMacro
using StaticArrays
using Artifacts 

#pkg used for intial conditions


#pkg used by the fluid evolution. 
using Interpolations
using DifferentialEquations
using Statistics
using JLD2


using LinearAlgebra
using LoopVectorization
using Surrogates
using StatsBase

using FastChebInterp


using TensorCast
using NumericalIntegration
using SimpleNonlinearSolve

using QuadGK


using HCubature
using UnPack

using PyPlot








const root_particle_lists=artifact"particle_lists"

#const root_kernel_particles=artifact"kernel"

const fmGeV= 1/0.1973261 
const invfmGeV= 1/fmGeV
const invfmGeV3=(1/fmGeV)^3
const fmGeV3=fmGeV^3






include("EquationOfState/thermodynamic.jl")
include("EquationOfState/EquationofStatetype.jl")
include("EquationOfState/TransportCoefficienttype.jl")


include("EquationOfState/IdealQCD.jl")
include("EquationOfState/FluiduMEoS.jl")

include("EquationOfState/charm_HRG.jl")
include("EquationOfState/simpletransport.jl")
include("EquationOfState/Analytic.jl")

include("EquationOfState/heavyquark.jl")


include("fluidevo/detector_struct.jl")
include("fluidevo/inverse_function.jl")



include("fluidevo/discretization.jl") 
include("fluidevo/picewisefunction.jl") 
include("fluidevo/discreteFields.jl")  

include("fluidevo/generating_kernels_3d.jl")
#include("fluidevo/generating_kernels.jl")
include("fluidevo/particle_dictionary.jl")


include("fluidevo/spectra_fastreso_dict_HQ.jl")
include("fluidevo/spectra.jl")

include("fluidevo/test_functions.jl")
include("fluidevo/map_profile.jl")
include("fluidevo/initial_fields.jl")

include("Matrix/1d_viscous_HQ_cilindrical_fugacity.jl")

include("wraps2.jl")
include("plotting_routines.jl")

const  detector_collection=(;ALICE=detector("ALICE",6.62,7.00,0.0757,"Pb_Pb"),
RHIC =detector("RHIC" ,7.,4.23,0.005968,"Au_Au"),
ALICE1 =detector("ALICE1",     6.62 ,    7.00 ,        0.0463,	"Pb_Pb")
)

const detector_dict=Dict(:ALICE=>detector(:ALICE,6.62,7.00,0.0757,:Pb_Pb),
:RHIC =>detector(:RHIC ,7.,4.23,0.005968,:Au_Au),
:ALICE1 =>detector(:ALICE1,     6.62 ,    7.00 ,        0.0463,	:Pb_Pb)
)


export detector_collection,detector_dict

export NDField, Fields, OriginInterval, CartesianDiscretization, DiscreteFileds
export set_array, set_array!, freeze_out_routine, jgemvavx!, oneshoot, SplineInterp, spectra_analitic
export spectra, spectra_lf, multiplicity, multiplicity_lf, DiscreteFields, TabulatedTrento, initialize_fields, FreezeOutResult
export get_profile, map_initial_profile, Profiles, Profiles2 #RunFluidum_hf, RunFluidum_array, save_to_h5, RunFluidum_lf, SetFluidProperties, FluidParameters



"""

    thermodynamic(T[,μ],x::EquationofState)
Give back the pressure, the pressure gradient and pressure hessian at given T and if specified at (T,μ) for given equation of state
"""
function thermodynamic
end

"""
    pressure(T[,μ],x::Equation of state)
Give back the pressure at given T and if specified at (T,μ) for given equation of state


"""
function pressure
end

"""
    pressure_derivative(T[,μ],Val(N)[,Val(M)],x::EquationofState)
Return the pressure gradient/hessian

# Arguments
- `T::Number`: Temperature
- `μ::Number`: Chemical potential
- `Val(N::Integer)`: N is the order of derivatives with respect to T ranging from 0 to 2
- `Val(M::Integer)`: M is the order of derivatives with respect to μ ranging from 0 to 2
- `x::EquationOfState`: equation of state

# Examples
Compute ∂_μ p for LatticeQCD at (T=0.5,μ=0.3)
```jldoctest
>julia LQCD=LatticeQCD()
LatticeQCD()
>julia pressure_derivative(0.5,0.3,Val(0),Val(1),LQCD)
0.024387044999999035
```
"""
function pressure_derivative
end

"""
    OneDPicewiseEquationOfState(eos1::EquationofState, int1::AbstractInterval,eos2::EquationOfState,int2::AbstractInterval)
Give back p_1/2(T), if T is in interval int1/2, else give back 0

"""
function OneDPicewiseEquationOfState
end

"""
    TwoDPicewiseEquationOfState(eos1::EquationofState, int1::AbstractInterval,eos2::EquationOfState,int2::AbstractInterval)
Give back p_1/2(T,μ), if (T,μ) is in interval int1/2, else give back 0
"""
function TwoDPicewiseEquationOfState
end

"""
    viscosity(T[,μ][,x::EquationOfState,th::Thermodynamic],y::TransportProperty)
Give back shear viscosity η/s and relaxation time τ_s 

# Arguments
- `T::Number`: Temperature
- `μ::Number`: Chemical potential
- `x::EquationOfState`: equation of state. Either this or thermodynamic has to be specified.
- `th::Thermodynamic`: thermodynamic. Either this or the equation of state has to be specified.
- `y::TransportProperty`: transport TransportProperty

# Examples
Compute the shear viscosity for T=0.3 GeV and μ=0.2 GeV for a hadron resonance gas
```jldoctest
>julia hrg=HadronResonaceGas()
HadronResonaceGas()
>julia vis=EquationsOfStates.SimpleShearViscosity(0.16,0.2)
SimpleShearViscosity{Float64}(0.16, 0.2)
>julia th=thermodynamic(0.3,0.2,hrg)
p(T,μ)=0.06585468437527325
∇p(T,μ)=(1.581024911648665, 0.07161301647333912)
∇²p(T,μ)=(30.990803012108447, 1.3978520817058184, 30.990803012108447)
>julia viscosity(0.3,0.2,hrg,vis)
0.0499163967709561
>julia viscosity(0.3,0.2,th,vis)
0.0499163967709561
```
"""
function viscosity
end

"""
    IdealQCD(Nf,Nc)
Set up fluid properties using the pressure of an ideal QCD gas with Nf flavours and Nc colors
"""
function IdealQCD
end

"""
    FluiduMEoS
Set up fluid properties using the pressure from a fit to HRG at low temperature and lattice QCD at high temperature
"""
function FluiduMEoS
end

"""
    HadronResonaceGas([,name_file,Maxmass,Minmass,condition])
Set up fluid properties using the pressure of a hadron resonance gas

# Arguments
- `name_file::String`: Name of file with particle information. Standart value is "EquationofState/particles.data"
- `Maxmass::Number`: Mass cut-off in GeV, heavier particles will be excluded. Standart value is 2.1 
- `Minmass::Number`: Mass cut-off in GeV, lighter particles will be excluded. Standart value is 0.05 
- `condition`: Custom condition to exclude particles
"""
function HadronResonaceGas
end


"""
    waleckacondition
Condition to exclude particles from hadron resonance gas already included in walecka model (ρ,ω,p^+,p^-)
"""
function waleckacondition
end

"""
    LatticeQCD
Set up fluid properties using the pressure from lattice QCD (paper ref)
"""
function LatticeQCD
end

"""
    WaleckaModel
Set up the fluid properties using the pressure from the waleckamodel
"""
function WaleckaModel
end

"""
    HadronResonaceGasNew
Set up the fluid properties from a hadron resonance gas based on the therminator particle file
"""
function HadronResonaceGasNew
end

"""
    SimpleBulkViscosity(ζs,Cζ)
Set up the bulk viscosity and corresponding relaxation time

# Arguments
- `ζs::Number`: Dimensionless value of maximum of bulk Bulk viscosity (ζ/s)_max, recommended value is 0.05
- `Cζ::Number`: Dimensionless value for magnitude of τ_bulk (∝ 1/C), recommended value is 15
"""
function SimpleBulkViscosity
end


"""
    SimpleShearViscosity(ηs,Cη)
Set up the shear viscosity and corresponding relaxation time

# Arguments
- `ηs::Number`: Value of shear viscosity over entropy ration η/s


TBW
"""
function SimpleShearViscosity
end

export thermodynamic,pressure, pressure_derivative , TwoDPicewiseEquationOfState ,OneDPicewiseEquationOfState, energy_density
export FluidProperties, viscosity, τ_shear,bulk_viscosity,τ_bulk,diffusion,τ_diffusion
export IdealQCD, FluiduMEoS, HadronResonaceGas,waleckacondition,LatticeQCD, HadronResonaceGasNew, HadronResonaceGas_ccbar
export SimpleBulkViscosity,SimpleShearViscosity,SimpleDiffusionCoefficient , ZeroViscosity, ZeroBulkViscosity
export Analytic, Gluing,Thermodynamic, EquationOfState, readresonancelist
export Heavy_Quark, HQdiffusion, free_charm, QGPViscosity,ZeroDiffusion

# type piracy 

Bessels.besselk0(d::Dual{T,V,N}) where {T,V,N} = Dual{T}(Bessels.besselk0(value(d)), -Bessels.besselk1(value(d)) * partials(d))

Bessels.besselk1(d::Dual{T,V,N}) where {T,V,N} = Dual{T}(Bessels.besselk1(value(d)),-( Bessels.besselk0(value(d))+Bessels.besselk1(value(d))/value(d)) * partials(d))

end
