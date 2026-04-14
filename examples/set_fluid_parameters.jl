#This code shows how to set the fluid parameters by using functions in heavy_quark.jl and in simple_transport.jl
using Fluidum
using Plots

#In order to inizialize the equation of state, we need the number of charm quarks created in the hard scattering. 
#In this example, we just use a plasuble number of the charm quarks for PbPb collisions at 5.02 TeV, with charm cross section computed from CTEQ
ccbar = 23.;
particle_list = string(Fluidum.root_particle_lists, "/OpenCharmParticleList_corrJS.txt");

#For the EoS, it is possible to specify the particle list and the number of charm quarks
eos = Heavy_Quark(particle_list, ccbar) 

#It is also possible to consider the default configuation, that takes the particle list specified in Fluidum.jl
eos = Heavy_Quark() 

#First example: without viscosity, bulk viscosity and diffusion
zero_viscosity = ZeroViscosity();
zero_bulk = ZeroBulkViscosity();
zero_diffusion = ZeroDiffusion();
zero_params = Fluidum.FluidProperties(eos,zero_viscosity,zero_bulk,zero_diffusion);

#second example: with viscosity, bulk viscosity and constant DsT. The values of the viscosity and bulk are inspired from literature arXiv:2308.16722v1[hep-ph]
#the value of the diffusion coefficient comes from arXiv:2302.08501 [hep-lat] evaluated at Tc, and maximized at the 
viscosity = QGPViscosity(0.1,0.2); 
bulk = SimpleBulkViscosity(0.1,15.0); 
constant_DsT = ConstDiffusion(0.24, 1.5);
params_const_Ds = Fluidum.FluidProperties(eos,viscosity,bulk,constant_DsT);

#third example: with viscosity, bulk viscosity and linear DsT. 
#The value of the slope and of the offset come from a linear fit in arXiv:2302.08501 [hep-lat], considering the central values 
slope = 1.765
offset = -0.159
linear_DsT = LinearDiffusion(slope, offset, 1.5);
params_linear_Ds = Fluidum.FluidProperties(eos,viscosity,bulk,linear_DsT);

