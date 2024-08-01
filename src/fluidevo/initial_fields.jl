


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

function initialize_fields(x::TabulatedTrento{A,B}, y::TabulatedTrento{C,D}, cent1::Integer, cent2::Integer, list; tau0 = 0.4, gridpoints=500,rmax=30,norm_temp=1, norm_coll=1, exp_tail = true) where {A,B,C,D}

    temperature_profile, nhard_profile = Profiles(x,y,cent1,cent2; radius = range(0, rmax, gridpoints), norm_temp = norm_temp, norm_coll = norm_coll, exp_tail = exp_tail)
    
    @show temperature_profile.(range(0,30,100))
    @show nhard_profile.(range(0,30,100)) 
    ccbar = quadgk(x->2*pi*x*tau0*nhard_profile(x),0,rmax,rtol=0.00001)[1]
    
    eos=HadronResonaceGas_ccbar(list, ccbar)

    fugacity(r) = log(nhard_profile(r)/(thermodynamic(temperature_profile(r),0.0,eos).pressure)) #.+ 0.0001
    @show fugacity.(range(0,30,100))
    @show ccbar
        
    ccbar_thermo,err= quadgk(x->2*pi*x*tau0*thermodynamic(temperature_profile(x),fugacity(x),eos).pressure,0,rmax,rtol=0.00001) #integrazione per ottenere il numero totale di charm
    @show ccbar_thermo

    oned_visc_hydro = Fluidum.HQ_viscous_1d()



    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    disc_fields = DiscreteFileds(oned_visc_hydro,disc,Float64) 
    phi=set_array((x)->temperature_profile(x),:temperature,disc_fields); #temperature initialization
    set_array!(phi,(x)->fugacity(x),:mu,disc_fields); #fugacity initialization
return DiscreteFields(disc,disc_fields,phi)
end

function initialize_fields(x::TabulatedTrento{A,B}, cent1::Integer, cent2::Integer, list; tau0 = 0.4, gridpoints=500,rmax=30,norm_temp=1, norm_coll=1, exp_tail = true) where {A,B,C,D}

    temperature_profile = Profiles(x,cent1,cent2; radius = range(0, rmax, gridpoints), norm_temp = norm_temp, norm_coll = norm_coll, exp_tail = exp_tail)
    
    @show temperature_profile.(range(0,30,100))
    
    oned_visc_hydro = Fluidum.viscous_1d()



    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    disc_fields = DiscreteFileds(oned_visc_hydro,disc,Float64) 
    phi=set_array((x)->temperature_profile(x),:temperature,disc_fields); #temperature initialization
    set_array!(phi,(x)->fugacity(x),:mu,disc_fields); #fugacity initialization
return DiscreteFields(disc,disc_fields,phi)
end







