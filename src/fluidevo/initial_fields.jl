struct DiscreteFields{A,B,C}    
    discretization::A
    discrete_field::B
    initial_field::C
end

@inline function HQ_viscous_1d() 
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:piphiphi),
    NDField((:even,),(:ghost,),:pietaeta),
    NDField((:even,),(:ghost,),:piB),
    NDField((:even,),(:ghost,),:mu), #fugacity
    NDField((:odd,),(:ghost,),:nur)
    )
end 

@inline function viscous_1d() 
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:piphiphi),
    NDField((:even,),(:ghost,),:pietaeta),
    NDField((:even,),(:ghost,),:piB)
    )
end 

function fug(T, nhard_profile, r, eos; rdrop = 12, m = 1.5)       
    if r<=rdrop
        return log(nhard_profile(r)/(thermodynamic(T(r),0.0,eos).pressure))
        else return log(nhard_profile(rdrop)/(thermodynamic(T(rdrop),0.0,eos).pressure))
    end            
end




function fug_fs(T, ncoll_profile, r, eos;dσ_QQdy, tau0, m = 1.5)       
    return log(density_polar(T, ncoll_profile,r;dσ_QQdy, tau0, m)/(thermodynamic(T,0.0,eos).pressure))            
end

function fug_fs_rdrop(T, ncoll_profile, r, eos;rdrop = 13, dσ_QQdy, tau0, m = 1.5)       
    if r <= rdrop 
        return fug_fs.(T(r), Ref(ncoll_profile), r, Ref(eos);dσ_QQdy,tau0,m=1.5)   
    else 
        ncoll_rdrop = x -> ncoll_profile(rdrop) 
        return fug_fs(T(rdrop), ncoll_rdrop, rdrop, eos;dσ_QQdy,tau0,m=1.5)
         
    end        
end



function nur_polar(T,ncoll_profile, fonll_profile,r; tau0, dσ_QQdy, m = 1.5) 
    nux_(x,y) = nux(T,ncoll_profile, fonll_profile,x,y; tau0, dσ_QQdy, m) 
    nur_polar = nux_(r,0) 
    return nur_polar
end


"""
Initialize the temperature, the fugacity, and the diffusion current 
"""
function initialize_fields(x::TabulatedData{A,B}, y::TabulatedData{C,D}, z::TabulatedData{E,F}, cent1::Integer, cent2::Integer, list; gridpoints=500, rmax=30, norm_temp=1, tau0 = 0.4, dσ_QQdy = 1, exp_tail = true,rdrop = 12, offset, m=1.5) where {A,B,C,D,E,F}
    temperature_profile, nhard_profile = Profiles(x,y,cent1,cent2; radius = range(0, rmax, gridpoints), norm_temp = norm_temp, norm_coll = 2/tau0*dσ_QQdy, exp_tail, offset)
    
    #@show temperature_profile.(range(0,rmax,gridpoints))
    #@show nhard_profile.(range(0,rmax,gridpoints))
    
    ccbar = quadgk(x->2*pi*x*tau0*nhard_profile(x),0,rmax,rtol=0.00001)[1]
    @show ccbar

    eos=HadronResonaceGas_ccbar(list, ccbar)
    fugacity_initial(r) = fug(temperature_profile, nhard_profile, r, eos;rdrop)   
        
    ccbar_thermo,err= quadgk(x->2*pi*x*tau0*thermodynamic(temperature_profile(x),fugacity_initial(x),eos).pressure,0,rmax,rtol=0.00001) #integrazione per ottenere il numero totale di charm
    @show ccbar_thermo
       
    fonll_profile = Profiles(z;norm_ = 1e-10, exp_tail = false)
    
    nur_profile(r) = nur_polar.(temperature_profile(r),Ref(nhard_profile), Ref(fonll_profile),r; tau0, dσ_QQdy, m)
    
    fugacity(r) = fug_fs_rdrop(temperature_profile, nhard_profile, r, eos;rdrop, dσ_QQdy, tau0, m = 1.5)  
        
    ccbar_thermo_fs,err= quadgk(x->2*pi*x*tau0*thermodynamic(temperature_profile(x),fugacity(x),eos).pressure,0,rmax,rtol=0.00001) #integration to obtain the total number of charm
    
    @show ccbar_thermo_fs

    oned_visc_hydro = Fluidum.HQ_viscous_1d()
    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    disc_fields = DiscreteFileds(oned_visc_hydro,disc,Float64) 
    phi=set_array((x)->temperature_profile(x),:temperature,disc_fields); #temperature initialization
    set_array!(phi,(x)->fugacity(x),:mu,disc_fields); #fugacity initialization
    
    set_array!(phi,(x)->nur_profile(x),:nur,disc_fields); #diffusion current initialization
return DiscreteFields(disc,disc_fields,phi)
end


"""
Use a step-like function for the collision profile if the collision profile is not given from trento
"""
function ncoll_step(r; Ncoll = 1, rdrop = 7, R = 7)       
    if r<=rdrop
        return  Ncoll/(pi*R^2)
        else return 0.01
    end            
end


"""
Initialize fields when only the temperature profile is given, using a step-like function for the collision profile
"""
function initialize_fields(x::TabulatedData{A,B},list; tau0 = 0.4, gridpoints=500,rmax=30,norm_temp=1, norm_coll=1, Ncoll = 1, rdrop = 7) where {A,B}
    temperature_profile = Profiles(x; norm_temp = norm_temp, radius = range(0, rmax, gridpoints)) 
    
    nhard_profile(r) = norm_coll*ncoll_step(r;Ncoll,rdrop)*(temperature_profile(r)/temperature_profile(0.0))^4 
    
    ccbar = quadgk(x->2*pi*x*tau0*nhard_profile(x),0,rmax,rtol=0.00001)[1]
    @show ccbar

    eos=HadronResonaceGas_ccbar(list, ccbar)
    fugacity(r) = fug(temperature_profile, nhard_profile, r, eos;rdrop)

    ccbar_thermo,err= quadgk(x->2*pi*x*tau0*thermodynamic(temperature_profile(x),fugacity(x),eos).pressure,0,rmax,rtol=0.00001) 
    @show ccbar_thermo
    
    n_therm(x) = thermodynamic(temperature_profile(x),fugacity(x),eos).pressure
    
    #@show n_therm.(range(0,rmax,gridpoints))
    
    #@show nhard_profile.(range(0,rmax,gridpoints))
    
    oned_visc_hydro = Fluidum.HQ_viscous_1d()
    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    disc_fields = DiscreteFileds(oned_visc_hydro,disc,Float64) 
    phi=set_array((x)->temperature_profile(x),:temperature,disc_fields); #temperature initialization
    set_array!(phi,(x)->fugacity(x),:mu,disc_fields); #fugacity initialization
    
return DiscreteFields(disc,disc_fields,phi)
end

