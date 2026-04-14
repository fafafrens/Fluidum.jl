struct DiscreteInitialFields{A,B,C}    
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
    NDField((:even,),(:ghost,),:α), 
    NDField((:odd,),(:ghost,),:nur)
    )
end

@inline function HQ_viscous_gamma_1d() 
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:piphiphi),
    NDField((:even,),(:ghost,),:pietaeta),
    NDField((:even,),(:ghost,),:piB),
    NDField((:even,),(:ghost,),:γ),
    NDField((:odd,),(:ghost,),:nur)
    )
end 

@inline function HQ_viscous_no_bulk_1d() 
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:piphiphi),
    NDField((:even,),(:ghost,),:pietaeta),
    NDField((:even,),(:ghost,),:α), 
    NDField((:odd,),(:ghost,),:nur)
    )
end 


@inline function HQ_viscous_BG_ideal_current_1d() 
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:piphiphi),
    NDField((:even,),(:ghost,),:pietaeta),
    NDField((:even,),(:ghost,),:piB),
    NDField((:even,),(:ghost,),:α)
    )
end 

@inline function HQ_viscous_BG_ideal_current_density_1d() 
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:piphiphi),
    NDField((:even,),(:ghost,),:pietaeta),
    NDField((:even,),(:ghost,),:piB),
    NDField((:even,),(:ghost,),:n)
    )
end 

@inline function HQ_ideal_1d() 
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:α),
    )
end 

@inline function HQ_ideal_br_p_density_1d() 
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:n),
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

function fug(T, nhard_profile, r, eos; rdrop = 20)       
    if r<=rdrop
        return log(nhard_profile(r)/(thermodynamic(T(r),0.0,eos).pressure))
        else return log(nhard_profile(rdrop)/(thermodynamic(T(rdrop),0.0,eos).pressure))
    end            
end

function gamma(T, nhard_profile, r, eos; rdrop = 20, m = 1.5)       
    if r<=rdrop
        @show nhard_profile(10)
        @show thermodynamic(T(10),0.0,eos).pressure
        return nhard_profile(r)/(thermodynamic(T(r),0.0,eos).pressure)
        else return nhard_profile(rdrop)/(thermodynamic(T(rdrop),0.0,eos).pressure)
    end            
end

function gamma_limit_one(T, nhard_profile, r, eos; rdrop = 20, m = 1.5)       
    if r<=rdrop
        return (nhard_profile(r)+r/30*(thermodynamic(T(r),0.0,eos).pressure-nhard_profile(r)))/(thermodynamic(T(r),0.0,eos).pressure)
        else return nhard_profile(rdrop)/(thermodynamic(T(rdrop),0.0,eos).pressure)
    end            
end


function fug_fs(T, ncoll_profile, r, eos;dσ_QQdy,tau_fs, tau0, m = 1.5)       
    return log(density_polar(T, ncoll_profile,r;dσ_QQdy, tau0,tau_fs, m)/(thermodynamic(T,0.0,eos).value))            
end

function fug_fs_rdrop(T, ncoll_profile, r, eos;rdrop = 13, dσ_QQdy, tau_fs, tau0, m = 1.5)       
    if r <= rdrop 
        return fug_fs.(T(r), Ref(ncoll_profile), r, Ref(eos);dσ_QQdy,tau0,tau_fs,m=1.5)   
    else 
        ncoll_rdrop = x -> ncoll_profile(rdrop) 
        return fug_fs(T(rdrop), ncoll_rdrop, rdrop, eos;dσ_QQdy,tau0,tau_fs,m=1.5)
         
    end        
end



"""
Initialize the temperature, the fugacity, and the diffusion current  
"""
function initialize_fields(x::TabulatedData{A,B}, y::TabulatedData{C,D}, z::TabulatedData{E,F}; cfg, m = 1.5) where {A,B,C,D,E,F}
    (; det, grid, tspan, tail_pars, cent) = cfg
    rmax = grid.rmax
    gridpoints = grid.gridpoints
    tau0 = tspan[1]
    
    temperature_profile, ncoll_profile = Profiles(x,y, cent; 
        radius = range(0, rmax, gridpoints),
        norm_x = cfg.temp_norm/tau0,
        norm_y = 2,
        tail_pars = tail_pars)

    fonll_profile = Profiles(z;norm = 1e-10, exp_tail = false)
    
    nhard_profile(r) = density_fs_total.(Ref(fonll_profile),Ref(ncoll_profile),r; tau0 = tau0, m = m) 
    
    ccbar_hard = quadgk(x->2*pi*x*tau0*nhard_profile(x),0,rmax,rtol=1e-5)[1]
    @show ccbar_hard

    Ncoll = nucl_radius(det)*quadgk(x->2*pi*x*nhard_profile(x)/(2/tau0*dσ_QQdy(det)),0,rmax,rtol=1e-5)[1]
    eos = Heavy_Quark(cfg.particle_list, ccbar_hard) 
    #fugacity(r) = fug_cut(temperature_profile, nhard_profile, r, eos.hadron_list; rmax = rmax)   
    
    fugacity(r) = fug(temperature_profile, nhard_profile, r, eos.hadron_list; rdrop = tail_pars.rdrop)   
    
    ccbar_thermo= quadgk(x->2*pi*x*tau0*thermodynamic(temperature_profile(x),fugacity(x),eos.hadron_list).pressure,0,rmax,rtol=1e-5)[1] #integration to obtain the total number of charm
    @show ccbar_thermo
        
    nur_profile(r) = nur.(temperature_profile(r),Ref(fonll_profile),Ref(ncoll_profile) ,r; tau0 = tau0, m)
        
    oned_visc_hydro = Fluidum.HQ_viscous_1d()
    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    disc_fields = DiscreteFields(oned_visc_hydro,disc,Float64) 
    phi=set_array((x)->temperature_profile(x),:temperature,disc_fields); #temperature initialization
    set_array!(phi,(x)->fugacity(x),:α,disc_fields); #fugacity initialization
    
    set_array!(phi,(x)->nur_profile(x),:nur,disc_fields); #diffusion current initialization
return DiscreteInitialFields(disc,disc_fields,phi), ccbar_thermo
end

"""
Initialize the temperature and the fugacity 
"""
function initialize_fields(x::TabulatedData{A,B}, y::TabulatedData{C,D}; cfg) where {A,B,C,D}
    (; det, grid, tspan, tail_pars, cent) = cfg
    rmax = grid.rmax
    gridpoints = grid.gridpoints
    tau0 = tspan[1]
    
    temperature_profile, nhard_profile = Profiles(x,y, cent; 
        radius = range(0, rmax, gridpoints),
        norm_x = cfg.temp_norm/tau0,
        norm_y = 2/tau0*dσ_QQdy(det),
        tail_pars = tail_pars)
    
    ccbar_hard = quadgk(x->2*pi*x*tau0*nhard_profile(x),0,rmax,rtol=1e-5)[1]  
    
    @show ccbar_hard
    eos = Heavy_Quark(cfg.particle_list, ccbar_hard)
    
    #fugacity(r) = fug_cut(temperature_profile, nhard_profile, r, eos.hadron_list; rmax = rmax)
    fugacity(r) = fug(temperature_profile, nhard_profile, r, eos.hadron_list; rdrop = tail_pars.rdrop)
     
    ccbar_thermo,err= quadgk(x->2*pi*x*tau0*thermodynamic(temperature_profile(x),fugacity(x),eos.hadron_list).pressure,0,rmax,rtol=1e-5) #integration to obtain the total number of charm
    @show ccbar_thermo
    
    oned_visc_hydro = Fluidum.HQ_viscous_1d()
    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    disc_fields = DiscreteFields(oned_visc_hydro,disc,Float64) 
    phi=set_array((x)->temperature_profile(x),:temperature,disc_fields); #temperature initialization
    set_array!(phi,(x)->fugacity(x),:α,disc_fields); #fugacity initialization
    
return DiscreteInitialFields(disc,disc_fields,phi), ccbar_hard
end


"""
Get initial temperature and fugacity profiles from tabulated data
"""
function get_initial_T_n_from_tabulated_data(x::TabulatedData{A,B}, y::TabulatedData{C,D}, cent1::Integer, cent2::Integer; grid_params, initial_params,dσ_QQdy, exp_tail = false) where {A,B,C,D}
    rmax = grid_params.rmax
    gridpoints = grid_params.gridpoints
    tau0 = initial_params.tau0

    temperature_profile, nhard_profile = Profiles(x,y,cent1,cent2; radius = range(0, rmax, gridpoints), norm_x = initial_params.norm/tau0,xmax_temp = 8, xmax_ncoll = 5, norm_y =2/tau0*dσ_QQdy, exp_tail)

    return Dict{Symbol,Function}(
        :temperature => (x -> temperature_profile(x)),
        :n           => (x -> nhard_profile(x)),
    )
end

"""
Get initial fugacity profile from temperature and density profiles
"""
function get_initial_α_from_T_and_n(temperature_profile, nhard_profile, eos, initial_params)

    fugacity = r -> fug(temperature_profile, nhard_profile, r, eos.hadron_list;initial_params.rdrop)  
    return Dict(
        :α => fugacity,
    )
end

function get_initial_γ_from_T_and_n(temperature_profile, nhard_profile, eos, initial_params)

    gamma = r -> gamma_limit_one(temperature_profile, nhard_profile, r, eos.hadron_list;initial_params.rdrop)  
    return Dict(
        :γ => gamma,
    )
end

"""
initialize fields given an dict of initialization functions
"""
function initialize_fields(init_fun_dict, field_initializer, grid_params)
    rmax = grid_params.rmax
    gridpoints = grid_params.gridpoints

    oned_visc_hydro = field_initializer
    disc = Fluidum.CartesianDiscretization(OriginInterval(gridpoints, rmax))
    disc_fields = Fluidum.DiscreteFields(oned_visc_hydro(), disc, Float64)

    phi = Fluidum.set_array((x) -> init_fun_dict[:temperature](x), :temperature, disc_fields) #temperature initialization

    try
        Fluidum.set_array!(phi, x -> init_fun_dict[:α](x), :α, disc_fields) #fugacity initialization
    catch
    end 
    try 
        Fluidum.set_array!(phi, x -> init_fun_dict[:γ](x), :γ, disc_fields) # gamma initialization
    catch
    end
    try
        Fluidum.set_array!(phi, x -> init_fun_dict[:n](x), :n, disc_fields) #density initialization
    catch
    end
    return Fluidum.DiscreteInitialFields(disc, disc_fields, phi)
end

"""
Initialize the temperature when only the temperature profile is given
"""
function initialize_fields_lf(x::TabulatedData{A,B}; cfg) where {A,B}
    (; grid, tspan, tail_pars, cent) = cfg
    rmax = grid.rmax
    gridpoints = grid.gridpoints
    tau0 = tspan[1]
    
    temperature_profile = Profiles(x,cent;
        radius = range(0, rmax, gridpoints), 
        norm = cfg.temp_norm/tau0, 
        tail_pars = tail_pars, 
        xmax = tail_pars.x_max_y,
        offset = tail_pars.offset_y,    
        temperature_flag = true)         

    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    oned_visc_hydro = Fluidum.viscous_1d()
    disc_fields = DiscreteFields(oned_visc_hydro,disc,Float64) 
    
    phi=set_array((x)->temperature_profile(x),:temperature,disc_fields); #temperature initialization
    
return DiscreteInitialFields(disc,disc_fields,phi)
end



"""
Use a step-like function for the collision profile if the collision profile is not given from trento
"""
function ncoll_step(r; Ncoll = 1653, rdrop, R)       
    if r<=rdrop
        return  Ncoll/(pi*R^2)
        else return 0.01
    end            
end


"""
Initialize fields when only the temperature profile is given, using a step-like function for the collision profile
"""
function initialize_fields(x::TabulatedData{A,B}; cfg) where {A,B}
    (; det, grid, tspan, tail_pars, cent) = cfg
    rmax = grid.rmax
    gridpoints = grid.gridpoints
    tau0 = tspan[1]
    
    temperature_profile = Profiles(x,cent;
        radius = range(0, rmax, gridpoints), 
        norm = cfg.temp_norm/tau0, 
        tail_pars = tail_pars,  
        xmax = tail_pars.x_max_y,
        offset = tail_pars.offset_y,    
        temperature_flag = true)         

    norm_y = 2/(tau0*σ_in(det))*dσ_QQdy(det)
    
    nhard_profile(r) = norm_y*ncoll_step(r; rdrop = tail_pars.rdrop, R = nucl_radius(det))*(temperature_profile(r)/temperature_profile(0.0))^4 
    ccbar_hard = quadgk(x->2*pi*x*tau0*nhard_profile(x),0,rmax,rtol=1e-5)[1]
    @show ccbar_hard

    eos = Heavy_Quark(cfg.particle_list, ccbar_hard)
    fugacity(r) = fug(temperature_profile, nhard_profile, r, eos.hadron_list;rdrop = tail_pars.rdrop)

    ccbar_thermo= quadgk(x->2*pi*x*tau0*thermodynamic(temperature_profile(x),fugacity(x),eos.hadron_list).pressure,0,rmax,rtol=1e-5)[1]
    @show ccbar_thermo
    
    n_therm(x) = thermodynamic(temperature_profile(x),fugacity(x),eos.hadron_list).pressure
    
    oned_visc_hydro = Fluidum.HQ_viscous_1d()
    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    disc_fields = DiscreteFields(oned_visc_hydro,disc,Float64) 
    phi=set_array((x)->temperature_profile(x),:temperature,disc_fields); #temperature initialization
    set_array!(phi,(x)->fugacity(x),:α,disc_fields); #fugacity initialization
    
return DiscreteInitialFields(disc,disc_fields,phi), ccbar_hard
end
