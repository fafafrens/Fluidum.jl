# If you have initial conditions tabulated data, you can use this script to test the causality of the evolution
using Fluidum
#define the grid
rmax = 20;
gridpoints=500;
grid_params = Fluidum.GridParameters(rmax, gridpoints);
grid_params

#centrality classes
cent1 = 0;
cent2 = 10;

#define the initial parameters
tspan = (0.4,10);
Norm = 90.; 
tau0 = tspan[1];
tau_fs = 0.1;
rdrop = 8.5;
initial_params = Fluidum.InitialParameters(Norm, tau0, tau_fs, rdrop);

#define the detector to use 
det = :ALICE;
det_name = detector_dict[det].name; 
dσ_QQdy = detector_dict[det].dσ_QQdy;
σ_in = detector_dict[det].σ_in;
entname = string("WeightFunction_", det_name, ".txt"); #path to entropy profile
ncollname = string("NcollProfile_", det_name, ".txt"); #path to ncoll profile
load_folder = string("/path/to/TabulatedData/");
#take the tabulated data 
ent = TabulatedData(load_folder * entname);
ncoll = TabulatedData(load_folder * ncollname);
fonll = TabulatedData(load_folder * fonllname);


disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
fields, ccbar = initialize_fields(ent, ncoll, cent1, cent2; grid_params, initial_params, dσ_QQdy);


eos = Heavy_Quark(readresonancelist(), ccbar)
viscosity = QGPViscosity(0.1,0.2); #or, ZeroViscosity();
bulk = SimpleBulkViscosity(0.083,15.0); #or, ZeroBulkViscosity();

diffusion_Ds = HQdiffusion(0.1, 1.5);
params_Ds=Fluidum.FluidProperties(eos,viscosity,bulk,diffusion_Ds);




field_results_Ds = Fluidum.oneshoot_debug(fields.discrete_field, Fluidum.matrxi1d_visc_HQ!, params_Ds, fields.initial_field, tspan; reltol = 1e-8);

# If you do NOT have initial conditions tabulated data, you can use this script to test the causality of the evolution
using Fluidum
#define equation of state and transport coefficients
ccbar = 30.
eos = Heavy_Quark(readresonancelist(), ccbar)

viscosity = QGPViscosity(0.1,0.2); #or, ZeroViscosity();
bulk = SimpleBulkViscosity(0.083,15.0); #or, ZeroBulkViscosity();

diffusion_Ds = HQdiffusion(0.1, 1.5);
params_Ds=Fluidum.FluidProperties(eos,viscosity,bulk,diffusion_Ds);
#define a radial grid
rmax = 20;
gridpoints=500;
grid_params = Fluidum.GridParameters(rmax, gridpoints);

disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
#define the name and number of fields, with appropriate parity
oned_visc_hydro = Fluidum.HQ_viscous_1d()
disc_fields = DiscreteFileds(oned_visc_hydro,disc,Float64) 
#define a example function for the temperature
function temperature(r)
       0.4*1/(exp(r/7)+1 )+0.0001
end

#set the temperature field
phi=set_array((x)->temperature(x),:temperature,disc_fields); 
set_array!(phi,(x)->-3. *temperature(x),:mu,disc_fields); 
    
tspan = (0.4,10);
field_results_Ds = Fluidum.oneshoot_debug(disc_fields, Fluidum.matrxi1d_visc_HQ!, params_Ds, phi, tspan; reltol = 1e-8);
