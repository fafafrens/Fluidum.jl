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

function fug(T, nhard_profile, r, eos; rdrop = 12)       
    if r<=rdrop
        return  log(nhard_profile(r)/(thermodynamic(T(r),0.0,eos).pressure))
        else return log(nhard_profile(rdrop)/(thermodynamic(T(rdrop),0.0,eos).pressure))
    end            
end

function fug_fs(T, ncoll_profile, x, y, eos; rdrop = 12,tau0, m = 1.5)       
    if sqrt(x^2+y^2)<=rdrop
        return log(density(ncoll_profile,x,y; tau0, m))/(thermodynamic(T(sqrt(x^2+y^2)),0.0,eos).pressure)
        else return log(density(ncoll_profile,rdrop/sqrt(2),rdrop/sqrt(2); tau0, m)/(thermodynamic(T(rdrop),0.0,eos).pressure))
    end            
end

function nux(T,ncoll_profile, fonll_profile,x,y; tau0, dσ_QQdy)
    ncoll_shift(px,py) = ncoll_fs(ncoll_profile,x,y,px,py;tau0=tau0) 
    fonll_interp(px,py) = fonll(fonll_profile,px,py; m = 1.5)
    eq_interp(px,py) = dσ_eq(T, sqrt(px^2+py^2);dσ_QQdy = dσ_QQdy)

    fonll_integral = hcubature(p->(p[1]*ncoll_shift(p[1],p[2])*fonll_interp(p[1],p[2])), rtol=0.001, [-20/sqrt(2), -20/sqrt(2)], [20/sqrt(2), 20/sqrt(2)])[1]
    eq_integral = hcubature(p->(p[1]*ncoll_shift(p[1],p[2])*eq_interp(p[1],p[2])), rtol=0.001, [-20/sqrt(2), -20/sqrt(2)], [20/sqrt(2), 20/sqrt(2)])[1]
    nux = 1/(2*dσ_QQdy)*(fonll_integral - eq_integral) 
    return nux
end

function nur_polar(T,ncoll_profile, fonll_profile,r; tau0 = 0.4, dσ_QQdy) 
    nux_(x,y) = nux(T,ncoll_profile, fonll_profile,x,y; tau0 = tau0, dσ_QQdy = dσ_QQdy) 
    nur_polar = nux_(r,0) 
    return nur_polar
end


"""
Initialize the temperature, the fugacity, and the diffusion current 
"""
function initialize_fields(x::TabulatedData{A,B}, y::TabulatedData{C,D}, z::TabulatedData{E,F}, cent1::Integer, cent2::Integer, list; gridpoints=500, rmax=30, norm_temp=1, tau0 = 0.4, dσ_QQdy = 1, exp_tail = true,rdrop = 12) where {A,B,C,D,E,F}
    temperature_profile, nhard_profile = Profiles(x,y,cent1,cent2; radius = range(0, rmax, gridpoints), norm_temp = norm_temp, norm_coll = 2/tau0*dσ_QQdy, exp_tail = exp_tail)
    
    #@show temperature_profile.(range(0,rmax,gridpoints))
    
    ccbar = quadgk(x->2*pi*x*tau0*nhard_profile(x),0,rmax,rtol=0.00001)[1]
    @show ccbar

    eos=HadronResonaceGas_ccbar(list, ccbar)
    fugacity(r) = fug(temperature_profile, nhard_profile, r, eos;rdrop)   
        
    ccbar_thermo,err= quadgk(x->2*pi*x*tau0*thermodynamic(temperature_profile(x),fugacity(x),eos).pressure,0,rmax,rtol=0.00001) #integrazione per ottenere il numero totale di charm
    @show ccbar_thermo
    
    fonll_profile = Profiles(z;norm_ = 1e-10, exp_tail = false)
    
    nur_profile(r) = nur_polar.(temperature_profile(r),Ref(nhard_profile), Ref(fonll_profile),r; tau0, dσ_QQdy)
    
    #nur_integral(r) = nur_int.(Ref(fonll_int), temperature_profile(r);dσ_QQdy = dσ_QQdy) 
    
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








