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
    NDField((:even,),(:ghost,),:mu), 
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



"""definition of fugacity"""
function fug(T, nhard_profile, r, eos; initial_params)       
    if r<=initial_params.rdrop
        return log(nhard_profile(r)/(thermodynamic(T(r),0.0,eos).pressure))
        else return log(nhard_profile(initial_params.rdrop)/(thermodynamic(T(initial_params.rdrop),0.0,eos).pressure))
    end            
end


"""
Initialize the temperature, the fugacity, and the diffusion current but rdrop  
"""
function initialize_fields(x::TabulatedData{A,B}, y::TabulatedData{C,D}, z::TabulatedData{E,F}, cent1::Integer, cent2::Integer; grid_params, initial_params, dσ_QQdy, exp_tail = true, m=1.5) where {A,B,C,D,E,F}
    rmax = grid_params.rmax
    gridpoints = grid_params.gridpoints
    tau0 = initial_params.tau0

    temperature_profile, ncoll_profile = Profiles(x,y,cent1,cent2; radius = range(0, rmax, gridpoints), norm_x = initial_params.norm/tau0, norm_y = 2, exp_tail)
    fonll_profile = Profiles(z;norm = 1e-10, exp_tail = false)
    
    nhard_profile(r) = density_fs_total.(Ref(fonll_profile),Ref(ncoll_profile),r; initial_params) 
    ccbar_hard = quadgk(x->2*pi*x*tau0*nhard_profile(x),0,rmax,rtol=0.00001)[1]
    ccbar_hard_after8 = quadgk(x->2*pi*x*tau0*nhard_profile(x),8,rmax,rtol=0.00001)[1]
    
    @show ccbar_hard
    @show ccbar_hard_after8

    Ncoll = 7*quadgk(x->2*pi*x*nhard_profile(x)/(2/tau0*dσ_QQdy),0,rmax,rtol=0.00001)[1]
    @show Ncoll
    eos = Heavy_Quark(readresonancelist(), ccbar_hard)
    fugacity(r) = fug(temperature_profile, nhard_profile, r, eos.hadron_list;initial_params)   
     
    
    ccbar_thermo,err= quadgk(x->2*pi*x*tau0*thermodynamic(temperature_profile(x),fugacity(x),eos.hadron_list).pressure,0,rmax,rtol=0.00001) #integration to obtain the total number of charm
    @show ccbar_thermo

    nur_profile(r) = nur.(temperature_profile(r),Ref(fonll_profile),Ref(ncoll_profile) ,r; initial_params, m)
    
    oned_visc_hydro = Fluidum.HQ_viscous_1d()
    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    disc_fields = DiscreteFileds(oned_visc_hydro,disc,Float64) 
    phi=set_array((x)->temperature_profile(x),:temperature,disc_fields); #temperature initialization
    set_array!(phi,(x)->fugacity(x),:mu,disc_fields); #fugacity initialization
    
    set_array!(phi,(x)->nur_profile(x),:nur,disc_fields); #diffusion current initialization
return DiscreteFields(disc,disc_fields,phi), ccbar_thermo
end



"""
Initialize the temperature and the fugacity but with rdrop 
"""
function initialize_fields(x::TabulatedData{A,B}, y::TabulatedData{C,D}, cent1::Integer, cent2::Integer; grid_params, initial_params, dσ_QQdy, exp_tail = true) where {A,B,C,D}
    rmax = grid_params.rmax
    gridpoints = grid_params.gridpoints
    tau0 = initial_params.tau0

    temperature_profile, nhard_profile = Profiles(x,y,cent1,cent2; radius = range(0, rmax, gridpoints), norm_x = initial_params.norm/tau0, norm_y = 2/tau0*dσ_QQdy, exp_tail)
        
    ccbar = quadgk(x->2*pi*x*tau0*nhard_profile(x),0,rmax,rtol=0.00001)[1]
    ccbar_out = quadgk(x->2*pi*x*tau0*nhard_profile(x),8,rmax,rtol=0.00001)[1]
    
    @show ccbar
    @show ccbar_out
    eos=Heavy_Quark(readresonancelist(), ccbar)
    fugacity(r) = fug(temperature_profile, nhard_profile, r, eos.hadron_list;initial_params)  
    ccbar_thermo,err= quadgk(x->2*pi*x*tau0*thermodynamic(temperature_profile(x),fugacity(x),eos.hadron_list).pressure,0,rmax,rtol=0.00001) #integration to obtain the total number of charm
    @show ccbar_thermo
        
    oned_visc_hydro = Fluidum.HQ_viscous_1d()
    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    disc_fields = DiscreteFileds(oned_visc_hydro,disc,Float64) 
    phi=set_array((x)->temperature_profile(x),:temperature,disc_fields); #temperature initialization
    set_array!(phi,(x)->fugacity(x),:mu,disc_fields); #fugacity initialization
    
return DiscreteFields(disc,disc_fields,phi), ccbar
end

"""
Initialize fields when only the temperature profile is given
"""
function initialize_fields(x::TabulatedData{A,B}, cent1::Integer, cent2::Integer; grid_params, initial_params, exp_tail = true) where {A,B}
    rmax = grid_params.rmax
    gridpoints = grid_params.gridpoints
    tau0 = initial_params.tau0
    
    temperature_profile = Profiles(x,cent1,cent2; radius = range(0, rmax, gridpoints), norm = initial_params.norm/tau0, exp_tail, temperature_flag = true)         


    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    oned_visc_hydro = Fluidum.viscous_1d()
    disc_fields = DiscreteFileds(oned_visc_hydro,disc,Float64) 
    
    phi=set_array((x)->temperature_profile(x),:temperature,disc_fields); #temperature initialization
    
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
    temperature_profile = Profiles(x; norm = norm_temp, radius = range(0, rmax, gridpoints)) 
    
    nhard_profile(r) = norm_coll*ncoll_step(r;Ncoll,rdrop)*(temperature_profile(r)/temperature_profile(0.0))^4 
    
    ccbar = quadgk(x->2*pi*x*tau0*nhard_profile(x),0,rmax,rtol=0.00001)[1]
    @show ccbar

    eos=HadronResonaceGas_ccbar(list, ccbar)
    fugacity(r) = fug(temperature_profile, nhard_profile, r, eos;rdrop)

    ccbar_thermo,err= quadgk(x->2*pi*x*tau0*thermodynamic(temperature_profile(x),fugacity(x),eos).pressure,0,rmax,rtol=0.00001) 
    @show ccbar_thermo
    
    n_therm(x) = thermodynamic(temperature_profile(x),fugacity(x),eos).pressure
    
    oned_visc_hydro = Fluidum.HQ_viscous_1d()
    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    disc_fields = DiscreteFileds(oned_visc_hydro,disc,Float64) 
    phi=set_array((x)->temperature_profile(x),:temperature,disc_fields); #temperature initialization
    set_array!(phi,(x)->fugacity(x),:mu,disc_fields); #fugacity initialization
    
return DiscreteFields(disc,disc_fields,phi)
end

"""
Initialize the temperature, the fugacity, and the diffusion current 
"""
function initialize_fields(x::TabulatedData{A,B}, y::TabulatedData{C,D}, cent1::Integer, cent2::Integer, list; gridpoints=500, rmax=30, norm_temp=1, tau0 = 0.4, dσ_QQdy = 1, exp_tail = true,rdrop = 12,offset = 0.01, m=1.5) where {A,B,C,D}
    temperature_profile, nhard_profile = Profiles(x,y,cent1,cent2; radius = range(0, rmax, gridpoints), norm_temp = norm_temp, norm_coll = 2/tau0*dσ_QQdy, exp_tail, offset=offset)
    
    ccbar = quadgk(x->2*pi*x*tau0*nhard_profile(x),0,rmax,rtol=0.00001)[1]
    @show ccbar

    eos=HadronResonaceGas_ccbar(list, ccbar)
    fugacity(r) = fug(temperature_profile, nhard_profile, r, eos;rdrop)   

    oned_visc_hydro = Fluidum.HQ_viscous_1d()
    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    disc_fields = DiscreteFileds(oned_visc_hydro,disc,Float64) 
    phi=set_array((x)->temperature_profile(x),:temperature,disc_fields); #temperature initialization
    set_array!(phi,(x)->fugacity(x),:mu,disc_fields); #fugacity initialization
    
return DiscreteFields(disc,disc_fields,phi)
end

