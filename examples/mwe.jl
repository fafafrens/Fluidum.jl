using Fluidum
using Plots
#define the grid
rmax = 20;
gridpoints=500;
grid_params = Fluidum.GridParameters(rmax, gridpoints);
grid_params
#define the detector to use 
det = :ALICE;
det_name = detector_dict[det].name; 
dσ_QQdy = detector_dict[det].dσ_QQdy;
σ_in = detector_dict[det].σ_in;

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


entname = string("WeightFunction_50k_PbPb502_", det_name, ".txt");
ncollname = string("NcollProfile_", det_name, ".txt");
fonllname = string("FONLL_", string(det), ".txt");
load_folder = string("/home/alice/Rossana_hq/TabulatedData/");
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

#fo_Ds = freeze_out_routine(fields.discrete_field, Fluidum.matrxi1d_visc_HQ!, params_Ds, fields.initial_field, tspan; Tfo = 0.156);
field_results_Ds = oneshoot(fields.discrete_field, Fluidum.matrxi1d_visc_HQ!, params_Ds, fields.initial_field, tspan; reltol = 1e-8);
iszero(Fluidum.oneshoot_debug(fields.discrete_field, Fluidum.matrxi1d_visc_HQ!, params_Ds, fields.initial_field, tspan; reltol = 1e-8))




n_func(t) = [thermodynamic(field_results_Ds(t)[1,i],field_results_Ds(t)[6,i],params_Ds.eos.hadron_list).pressure for i in eachindex(field_results_Ds(0.4)[1,:])]
plot(log.(field_results_Ds(0.4)[1,:]))
plot!(field_results_Ds(1.4)[7,:])
plot!(field_results_Ds(4.9)[7,:])
#plot!(field_results_Ds(5.4)[7,:])
#plot!(field_results_Ds(8.4)[7,:])
n_func(0.4)
plot(n_func(0.4))
plot!(n_func(3.4))

plot(field_results_Ds(0.4)[6,:])
plot!(field_results_Ds(1.4)[6,:])
plot!(field_results_Ds(4.9)[6,:])
#plot!(field_results_Ds(8.4)[6,:])


r, entropy_profile = get_profile(ent, cent1, cent2; norm = Norm)
    
temperature_profile = Fluidum.InverseFunction(x->pressure_derivative(x,Val(1),FluiduMEoS())).(entropy_profile)  
temperature_funct = Fluidum.linear_interpolation(r, temperature_profile; extrapolation_bc=Fluidum.Flat()) #from 9.9 fm it's flat
    
temp_exp = Fluidum.exponential_tail_pointlike.(Ref(temperature_funct), range(0, rmax, gridpoints); xmax = 8, offset = 0.01)
   # temp_exp = exponential_tail_pointlike.(Ref(temperature_funct), radius; xmax = 8, offset = 0.01)
temp_exp_funct = Fluidum.linear_interpolation(range(0, rmax, gridpoints), temp_exp; extrapolation_bc=Fluidum.Flat()) 

temperature_funct_grid = [temperature_funct(i) for i in range(0, rmax, gridpoints)]
temperature_exp_funct_grid = [temp_exp_funct(i) for i in range(0, rmax, gridpoints)]

plot(range(0, rmax, gridpoints),log.(temperature_funct_grid))
plot!(range(0, rmax, gridpoints),log.(temp_exp))
plot!(range(0, rmax, gridpoints),log.(temperature_exp_funct_grid))

#ncoll profile
r, ncoll_profile = get_profile(ncoll, cent1, cent2; norm = 1)  
ncoll_funct = Fluidum.linear_interpolation(r, ncoll_profile; extrapolation_bc=Fluidum.Flat()) #from 9.9 fm it's flat
   
ncoll_exp = Fluidum.exponential_tail_pointlike.(Ref(ncoll_funct), range(0, rmax, gridpoints);xmax = 4.5, offset=0.005)
        
ncoll_exp_funct = Fluidum.linear_interpolation(range(0, rmax, gridpoints), ncoll_exp; extrapolation_bc=Fluidum.Flat()) 
ncoll_funct_grid = [ncoll_funct(i) for i in range(0, rmax, gridpoints)]
ncoll_exp_funct_grid = [ncoll_exp_funct(i) for i in range(0, rmax, gridpoints)]
  
plot(range(0, rmax, gridpoints),log.(ncoll_funct_grid))
plot!(range(0, rmax, gridpoints),log.(ncoll_exp))
plot!(range(0, rmax, gridpoints),log.(ncoll_exp_funct_grid))
