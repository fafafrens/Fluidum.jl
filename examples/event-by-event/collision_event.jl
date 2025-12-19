
using Fluidum
using MonteCarloGlauber
include("MCglauber.jl")
include("observables.jl")
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
gridpoints=200
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
    mult, x_com, y_com = center_of_mass(event,Nr=500, Nth=500)
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

struct ObservableResult
    GlauberMultiplicity::Float64
    pTlist::Vector{Float64}
    vn::Array{Float64,4} # (cos, sin, denom) x length(pTlist) x length(wavenum_m) x length(species_list)
end

notempty = run_event(
    participants,
    twod_visc_hydro_discrete,
    100;
    pTlist = collect(0.1:0.2:2.0)
)

notempty.pTlist
notempty.vn
notempty.GlauberMultiplicity

empty = run_event(
    participants,
    twod_visc_hydro_discrete,
    0.1;
    pTlist = collect(0.1:0.2:2.0)
)

empty.pTlist
empty.vn
empty.GlauberMultiplicity

length(empty.points)
length(notempty.points)
using HDF5


result_file = pwd()*"/examples/event-by-event/results_auto.h5"

function run_event_by_event(Nev,namefile)
    h5open(namefile, "w") do file
        for i in 1:Nev
            result = run_event(
                participants,
                twod_visc_hydro_discrete,
                norm;
                pTlist = collect(0.5:0.5:4.0)
            )

           
            dset_name = "/run_$i"
            write(file, dset_name, result)
            println(i)
           
        end
    end
    close(file)
end

run_event_by_event(10,result_file)

result_array = []
function run_event_by_event(result_array,Nev)
        for i in 1:Nev
            result = run_event(
                participants,
                twod_visc_hydro_discrete,
                norm;
                pTlist = collect(0.5:0.5:4.0)
            )
            push!(result_array,result)
            println(i)
        end
    end


example = run_event(
    participants,
    twod_visc_hydro_discrete,
    norm;
    pTlist = collect(0.1:0.2:2.0)
)


run_event_by_event(10)
result_array[1].vn
pTlist = result_array[1].pTlist



M_species(result_array[10],[pion,D0])
g_species_ptbin(result_array[10],[pion,D0])[1,2]

function qvec(result_single_event,species_list,wavenum_list)
    qvec_result = zeros(length(wavenum_list))
    pTlist = result_single_event.pTlist
    vn = result_single_event.vn

    for wavenum in eachindex(wavenum_list)
        for k in eachindex(species_list)
            for pt_idx in eachindex(pTlist)
                qvec_result[wavenum]+=sqrt(vn[1,pt_idx,wavenum,k]^2+vn[2,pt_idx,wavenum,k]^2)/vn[3,pt_idx,wavenum,k]*g_species_ptbin(result_single_event,species_list)[pt_idx,k]
            end
        end
    end
    return qvec_result
end

qvec(result_array[10],[pion,D0],[2,3])


function v_wavenum(result_array,species_list,wavenum_list)
    
    pTlist = result_array[1].pTlist
    vm_final= zeros(length(pTlist),length(wavenum_list),length(species_list))

    for pt_idx in eachindex(pTlist)
    for wavenum in eachindex(wavenum_list)
    for k in eachindex(species_list)
        for i in eachindex(result_array)

       
        vn = result_array[i].vn

        q_vector = qvec(result_array[i],species_list,wavenum_list)

        vm_final[pt_idx,wavenum,k] += 1/Nev*sqrt(vn[1,pt_idx,wavenum,k]^2+vn[2,pt_idx,wavenum,k]^2)/vn[3,pt_idx,wavenum,k]*q_vector[wavenum]
        

        end
    end
    end
    end
    return vm_final
end

vns = v_wavenum(result_array,[pion,D0],[2,3])

scatter(pTlist,vns[:,1,1],label="v2 pion")
scatter!(pTlist,vns[:,1,2],label="v2 D0")
scatter!(pTlist,vns[:,2,1],label="v3 pion")
scatter!(pTlist,vns[:,2,2],label="v3 D0")

Nev = length(result_array)

for pt_idx in eachindex(pTlist)

    for i in eachindex(result_array)
        qvec_pion[pt_idx] = sqrt(result_array[i].vn[1,pt_idx,1,1]^2+result_array[i].vn[2,pt_idx,1,1]^2)/result_array[i].vn[3,pt_idx,1,1]
        qvec_D0[pt_idx] =sqrt(result_array[i].vn[1,pt_idx,1,2]^2+result_array[i].vn[2,pt_idx,1,2]^2)/result_array[i].vn[3,pt_idx,1,2]
       
        v2_pion[pt_idx]=v2_pion[pt_idx]+1/Nev*qvec_pion[pt_idx]*qvec
        v2_D0[pt_idx]=v2_D0[pt_idx]+1/Nev*qvec_D0[pt_idx]*qvec
    end
end

v2_pion

using Plots
scatter(pTlist,v2_pion)
scatter!(pTlist,v2_D0)

file = h5open("results.h5", "r")
vn = Array{Float64}(undef,1,3,10,2)

pTlist = collect(0.1:0.2:2.0)
elliptic_flow = zeros(length(pTlist))
triangular_flow = zeros(length(pTlist))
for i in 1:17
    dset_name = "/run_$i"
    data = read(file, dset_name)
    vn[i,:,:,:] = data
end
close(file)



    for (i,pT) in enumerate(pTlist)
        sum_cos = 0.
        sum_sin = 0.
        sum_den = 0.
            
        for j in 1:17
            sum_cos = vn[j,1,i,1]
            sum_sin = vn[j,2,i,1]
            sum_den = vn[j,3,i,1]
            elliptic_flow[i] += sqrt(sum_cos^2 + sum_sin^2)/sum_den
        end
    end


elliptic_flow



    for (i,pT) in enumerate(pTlist)
            sum_cos += vn[j,1,i,2]
            sum_sin += vn[j,2,i,2]
            sum_den += vn[j,3,i,2]
            for j in 1:17
            sum_cos = 0.
            sum_sin = 0.
            sum_den = 0.
            triangular_flow[i] = (sum_cos^2 + sum_sin^2)/sum_den
        end
    end

    using Plots
scatter(pTlist, elliptic_flow, label="v2")

scatter(pTlist, triangular_flow, label="v3")