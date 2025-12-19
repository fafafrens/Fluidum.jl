"""
nhard_profile(x,y,ncoll;rmax=8)
distributes ncoll uniformly in a circle of radius rmax,0 outside
"""
#missing the tau0 factor, 2 for charm + anticharm, sigma in and dsigma dy for this collision energy
function nhard_profile(x,y,ncoll;rmax=8)
  #  if x^2 + y^2 < rmax^2
        #return ncoll/(pi*rmax^2)
        return exp(- (x^2 + y^2)/(2*(rmax/2)^2))*ncoll/(2*pi*(rmax/2)^2)
        #return 1/(1+exp(x^2+y^2))*ncoll + 0.01
  #  else
  #      return 0.00
  #  end
end

function fug_(profile, ncoll, eos, discretization)   
    fugacity(x,y,T)=log(nhard_profile(x, y, ncoll)/(thermodynamic(T,0.0,eos.hadron_list).pressure))
    fug_map = map(zip(discretization.grid,profile)) do x
            return fugacity(x[1][1],x[1][2],x[2])

    end
    #fug_exp = [Fluidum.exponential_tail_pointlike(fug_func, x, y; offset=0.0015) for x in xgrid, y in ygrid]
    #fug_exp_func = linear_interpolation((xgrid,ygrid), fug_exp; extrapolation_bc=Flat()) 
   
    return fug_map
end

function trento_event_discrete(profile;norm=10, xgrid = -10:0.5:10, ygrid = -10:0.5:10)
    
    entropy_profile_ext = linear_interpolation((xgrid,ygrid),norm.*profile,extrapolation_bc=Flat())
    return map(Iterators.product(xgrid,ygrid)) do I
    entropy_profile_ext(I...)
    end
end



function trento_event_eos(profile;norm=10,exp_tail = false, xgrid = -10:0.5:10, ygrid = -10:0.5:10, kwarg...)
    
   # entropy_profile_ext = linear_interpolation((xgrid,ygrid),norm.*profile,extrapolation_bc=Interpolations.Flat())
   # entropy_profile = map(Iterators.product(xgrid,ygrid)) do I
   # entropy_profile_ext(I...)
   # end

    temperature_profile = InverseFunction(x->pressure_derivative(x,Val(1),FluiduMEoS())).(norm.*profile)
   # temperature_funct = linear_interpolation((xgrid,ygrid), temperature_profile; extrapolation_bc=Flat())
    
   if exp_tail==true
   
    temp_exp = [Fluidum.exponential_tail_pointlike(temperature_funct, x, y; offset=0.01) for x in xgrid, y in ygrid]
    temp_exp_funct = linear_interpolation((xgrid,ygrid), temp_exp; extrapolation_bc=Flat()) 
        return  temp_exp_funct
    else
        return temperature_profile
    end
  
end

function runGlauber(Nevents,centrality,n1,n2,w,s_NN,k,p)

participants=Participants(n1,n2,w,s_NN,k,p);
event=rand(participants,Nevents);

batches, CoM=MonteCarloGlauber.centralities_selection_CoM(event,centrality)

batches[1];
end

function set_array(phi::AbstractArray,express::S,disc::DiscreteFileds{T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,S,N_field2}) where {S,T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,N_field2}
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

function set_array!(array,phi::AbstractArray,express::S,disc::DiscreteFileds{T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,S,N_field2}) where {S,T,total_dimensions,space_dimension,N_field,Sizes_ghosted,Sizes,Lengths,DXS,M,N_field2}
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
