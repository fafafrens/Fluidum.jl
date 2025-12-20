
using Fluidum
using MonteCarloGlauber
using MuladdMacro
using OhMyThreads
include("MCglauber.jl")
include("observables.jl")
include("hdf5_io.jl")
# the convention here are T, ux, uy, piyy, pizz, pixy, piB this has to match with the matrix defined
twod_visc_hydro=Fields(
NDField((:even,:ghost),(:even,:ghost),:temperature),
NDField((:ghost,:ghost),(:ghost,:ghost),:ux),
NDField((:ghost,:ghost),(:ghost,:ghost),:uy),
NDField((:ghost,:ghost),(:even,:ghost),:piyy),
NDField((:even,:ghost),(:even,:ghost),:pizz),
NDField((:ghost,:ghost),(:ghost,:ghost),:pixy),
NDField((:even,:ghost),(:even,:ghost),:piB),
NDField((:even,:ghost),(:even,:ghost),:mu),
NDField((:ghost,:ghost),(:ghost,:ghost),:nux),
NDField((:ghost,:ghost),(:ghost,:ghost),:nuy)
)

#we define a 2 cartesian grid form -25 to 25 50 point each dimension 
gridpoints=100
xmax = 25.
discretization=CartesianDiscretization(Fluidum.SymmetricInterval(gridpoints,xmax),Fluidum.SymmetricInterval(gridpoints,xmax))
# we prepare the field with the discretization
twod_visc_hydro_discrete=DiscreteFileds(twod_visc_hydro,discretization,Float64)

#nuclear parameters
n1= TabulatedEvent(pwd()*"/examples/event-by-event/NLEFT_dmin_0.5fm_positiveweights_O.h5")
n2= n1
w= 1
s_NN=5000
k=1
p=0.

#entropy normalization
norm = 90
participants=Participants(n1,n2,w,s_NN,k,p)

pion = particle_simple(0.13957,1,0)
D0 = particle_simple(1.86483,1,1)

function run_event(participants,twod_visc_hydro_discrete,norm;pTlist=collect(0.1:0.5:2.0),eta_p=0.0,
    wavenum_m=[2,3],species_list = [pion,D0])
    discretization=twod_visc_hydro_discrete.discretization
    #create event
    event=rand(participants);
    #compute center of mass
    mult, x_com, y_com = center_of_mass(event,100,50)
    xcm= x_com/mult
    ycm= y_com/mult   
    profile=map(discretization.grid) do y
            y = y.+(xcm,ycm)
            event(y...)
    end
    ncoll_event=event.n_coll

    #set up fluid properties
    dσ_QQdy = 0.05688 #in mb FONLL
    ccbar = ncoll_event*dσ_QQdy/70/0.4
    eos = Heavy_Quark(readresonancelist(), ccbar)
    viscosity=QGPViscosity(0.2,0.2)
    bulk=SimpleBulkViscosity(0.05,15.0)
    diffusion=HQdiffusion(0.1,1.5)

    fluidproperty=FluidProperties(eos,viscosity,bulk,diffusion)
    
    #setup fields
    temperature_func = trento_event_eos(profile,norm=norm,exp_tail=false)
    fug_func = fug_(temperature_func,ncoll_event, eos, discretization)   
    phi=set_array(temperature_func.+0.01,:temperature,twod_visc_hydro_discrete);
    set_array!(phi,fug_func,:mu,twod_visc_hydro_discrete);

    tspan=(0.4,30.)
    Tfo=0.156
    result=Fluidum.isosurface(twod_visc_hydro_discrete,Fluidum.matrxi2d_visc_HQ!,fluidproperty,phi,tspan,:temperature,Tfo)
    cha=Fluidum.Chart(Fluidum.Surface(result[:surface]),(t,x,y)->Fluidum.SVector{2}(atan(t,hypot(y,x)),atan(y,x)))
    if length(cha.points)==0
        return ObservableResult(mult,zeros(length(pTlist)),Array{Float64}(undef,3,length(pTlist),length(wavenum_m),length(species_list)))
    end
    fo_bg=Fluidum.freezeout_interpolation(cha,sort_index=2,ndim_tuple=50)

    #run observables
    vn = dvn_dp_list_delta(fo_bg,species_list, pTlist, eta_p, wavenum_m; eta_min=-5.0, eta_max=5.0)
   # spectra = dn_dpTdetap(fo_bg,m,pTlist, eta_p; eta_min=-5.0, eta_max=5.0)
   # writedlm("dvn_dp_list_result.txt", vn)
return ObservableResult(mult,pTlist,vn.u)
end


function run_event_by_event(Nev)
        tmap(1:Nev) do i
            println("Running event $i / $Nev")
            result = run_event(
                participants,
                twod_visc_hydro_discrete,
                norm;
                pTlist = collect(0.5:0.5:4.0)
            )
            result
        end
 end



Nev = 100_000
checkpoint_interval = 1000
checkpoint_file = "event_by_event_results.h5"

for batch in 1:ceil(Int, Nev / checkpoint_interval)
    batch_size = min(checkpoint_interval, Nev - (batch - 1) * checkpoint_interval)
    println("Running batch $batch: $batch_size events")
    data = run_event_by_event(batch_size)
    append_to_h5(checkpoint_file, data)
    println("Checkpoint saved. Total events: $(batch * checkpoint_interval)")
end


