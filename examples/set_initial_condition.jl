#This code shows how to set the initial fields by using functions written in initial_fields.jl and map_profile.jl
using Fluidum
using Plots

"""We define here some plausible analytic profiles for entropy, density of binary collision and charm production cross section
In practice these profiles should come from Trento (entropy and density of binary collisions) and pQCD fonll calculations (charm production cross section)."""
function entropy(r)
       2/(exp((r-5)/0.75)+1)
end

function ncoll(r)
    4.2/(exp((r-4)/0.8)+1)
end


function fonll(pt)
       1.7e7*(exp(-(pt/2.5)^2))
end

#In this case, we need to define the radial and momentum range. In normal applications, these are already given from Trento and fonll calculations.
rmax = 30.
gridpoints = 300

pt_max = 20.
pt_points = 100

r = collect(range(0,rmax,gridpoints))
pt = collect(range(0.1,pt_max,pt_points))

"""convert the analytic profiles into the format TabulatedData, which is the format used for the initializations. 
Repeat the entropy and binary collision profiles in 100 columns, to mimic the output of trento with different centrality classes """
tab_entropy = TabulatedData(r, repeat(reshape(entropy.(r), :, 1), 1, 100))
tab_ncoll = TabulatedData(r, repeat(reshape(ncoll.(r), :, 1), 1, 100))
tab_fonll = TabulatedData(pt, reshape(fonll.(pt), :, 1))


"""initial configuration parameters. 
Available detector names: ALICE and RHIC, at energies 5.02 TeV and 0.200 TeV, respectively.
Available colliding nuclei: PbPb and AuAu.
Available PDFs: CTEQ, NNPDF and experimental fits, considering shadowing effects.

norm refers to the temperature normalization, tau0 is the QGP formation time and rdrop the radius from which we consider the fugacity to be constant, when we apply the exponential tail 

The offset in x is related to the temperature, while the one in y is related to the ncoll profile. 
The xmax parameters determine where the exponential tail starts to be applied.  
"""
cfg = RunConfig(
    det = Detector(;det_name = "ALICE", energy = 5.02, nuclei = "PbPb", PDF = "CTEQ"), #detector specifications  
    grid = GridParameters(rmax = rmax, gridpoints = gridpoints), #grid parameters
    temp_norm = 90.,
    tail_pars = ExpTail(is_tail = true, rdrop = 8., offset_x = 0.005, offset_y = 0., x_max_x = 8., x_max_y = 5.), #exponential tail parameters
    cent = Centrality(0,10), #centrality classes
    tspan = (0.4,12.), 
    Tfo = 0.156, #freeze-out temperature
    particle_list = string(Fluidum.root_particle_lists, "/OpenCharmParticleList_corrJS.txt"), 
)


"""Otherwise, it is possible to define a RunConfig with default PDF, grid, exponential tail and particle list.
The default values are reported in initial_parameters.jl"""
cfg2 = RunConfig(
    det = Detector(det_name = "ALICE", energy = 5.02, nuclei = "PbPb"),
    temp_norm = 90.,
    tspan = (0.4, 12.),
    cent = Centrality(0,10),
    Tfo = 0.156,
)


#initialize the temperature profile only
fields = Fluidum.initialize_fields_lf(tab_entropy; cfg = cfg);

#initialize the temperature profile and the HQ fugacity, using a step-like function for the collision profile
fields, ccbar = Fluidum.initialize_fields(tab_entropy; cfg = cfg);

#initialize the temperature profile and the HQ fugacity 
fields, ccbar = Fluidum.initialize_fields(tab_entropy, tab_ncoll; cfg = cfg);

#initialize the temperature profile, the HQ fugacity and the initial diffusion current
fields, ccbar = Fluidum.initialize_fields(tab_entropy, tab_ncoll, tab_fonll; cfg = cfg);

#get the initial fields and the discretization from the initialization
phi = fields.initial_field;
discretization = fields.discretization;

#get the initial grid
x = [discretization.grid[i][1] for i in 2:lastindex(discretization.grid)-1];
 
#evaluate the temperature, the fugacity and the diffusion current (the three fields that were initialized)
T(x) = phi[1,eachindex(x)]
fug(x) = phi[6,eachindex(x)]
nu(x) = phi[7,eachindex(x)]


#plot the resulting profiles
plot(x, T(x), xlabel = "r [fm]", ylabel = "T [GeV]")
plot(x, fug(x), xlabel = "r [fm]", ylabel = "fugacity")
plot(x, nu(x), xlabel = "r [fm]", ylabel = "diffusion current [fm^-3]")
