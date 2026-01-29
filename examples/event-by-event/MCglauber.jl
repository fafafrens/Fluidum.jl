"""
nhard_profile(x,y,ncoll;rmax=8)
distributes ncoll with a gaussian
"""
function nhard_profile(x,y,ncoll;rmax=8,norm = 2*0.4087/(70.)/0.4)
        norm*exp(- (x^2 + y^2)/(2*(rmax/2)^2))*ncoll/(2*pi*(rmax/2)^2) #2d-gaussian with σ=rmax/2      
end

function charm_number_hard(ncoll;xmax=20.,norm = 2*0.4087/(70.)/0.4)
    domain = ([-xmax,-xmax],[xmax,xmax])
    function f(u,p) 
        ncoll = p
        nhard_profile(u[1],u[2],ncoll,norm=norm)
    end
    p = (ncoll)
    prob = IntegralProblem(f,domain,p)
    result = solve(prob, HCubatureJL(), reltol=1e-3, abstol=1e-6)
    return result
end


function fug_(profile, ncoll, eos, discretization; norm = 2*0.4087/(70.)/0.4)   
    fugacity(x,y,T)=log(nhard_profile(x, y, ncoll; norm = norm)/(thermodynamic(T,0.0,eos.hadron_list).pressure))
    fug_map = map(zip(discretization.grid,profile)) do x
            return fugacity(x[1][1],x[1][2],x[2])

    end
   
    return fug_map
end

function trento_event_discrete(profile;norm=10, xgrid = -10:0.5:10, ygrid = -10:0.5:10)
    
    entropy_profile_ext = linear_interpolation((xgrid,ygrid),norm.*profile,extrapolation_bc=Flat())
    return map(Iterators.product(xgrid,ygrid)) do I
    entropy_profile_ext(I...)
    end
end



function trento_event_eos(profile;norm=10,exp_tail = false, xgrid = -10:0.5:10, ygrid = -10:0.5:10, kwarg...)
 

    temperature_profile = InverseFunction(x->pressure_derivative(x,Val(1),FluiduMEoS())).(norm.*profile)
    
   if exp_tail==true
   
    temp_exp = [Fluidum.exponential_tail_pointlike(temperature_funct, x, y; offset=0.01) for x in xgrid, y in ygrid]
    temp_exp_funct = linear_interpolation((xgrid,ygrid), temp_exp; extrapolation_bc=Flat()) 
        return  temp_exp_funct
    else
        return temperature_profile
    end
  
end

function set_array(phi::AbstractArray,express::S,disc::DiscreteFields{T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,S,N_field2}) where {S,T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,N_field2}
    array=Fluidum.get_array(disc)
    i=Fluidum.get_index(express,disc.fields)
    for I in disc.index_structure.interior
        #point=disc.discretization[I]
        point=disc.discretization.grid[I]
        array[i,I]= phi[I]
    end
    Fluidum.boundary_condition!(array,disc)
    return array
end

function set_array!(array,phi::AbstractArray,express::S,disc::DiscreteFields{T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,S,N_field2}) where {S,T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,N_field2}
    #array=zeros(T,(N_field,Sizes_ghosted...))
    i=Fluidum.get_index(express,disc.fields)
    for I in disc.index_structure.interior
        #point=disc.discretization[I]
        point=disc.discretization.grid[I]
        array[i,I]= phi[I]
    end

    Fluidum.boundary_condition!(array,disc)

    return array
end
