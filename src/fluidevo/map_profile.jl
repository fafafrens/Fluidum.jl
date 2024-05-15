


abstract type AbstractInitialCondition end

struct TabulatedTrento{A<:AbstractArray, B<:AbstractArray} <: AbstractInitialCondition
    radius::A
    profile::B
end

function TabulatedTrento(fn::A) where {A<:AbstractString}
    file = readdlm(fn)
    TabulatedTrento(file[:,1], file[:,2:end])
end   


function centrality_bin(x::TabulatedTrento{A,B}) where {A,B}
    n_columns = size(x.profile)[2]
    bin_size = div(100, n_columns)
end

function get_profile(x::TabulatedTrento{A,B}, cent1::Integer, cent2::Integer; norm = 1) where {A,B}
    r = x.radius
    delta_bin = centrality_bin(x) 
    idx1 = div(cent1,delta_bin)+1
    idx2 = div(cent2,delta_bin)
    mean_profile = @views dropdims(Statistics.mean(x.profile[:,idx1:idx2]; dims = 2); dims=2)
    
    if mod(cent1, delta_bin) != 0 || mod(cent2, delta_bin) != 0 
        @warn string("centrality chosen is not a multiple of the bin size. Using ", (idx1-1)*delta_bin, "-", idx2*delta_bin, " instead")       
    end 
    
    return r, norm.*mean_profile
end


function exponential_tail_pointlike(function_profile,x,xmax; offset = 0.015)
    a = ForwardDiff.derivative(function_profile, xmax) / function_profile(xmax)
    b = log(function_profile(xmax))
    if x-xmax < 0
        return function_profile(x)
    else  
        return exp(a*(x-xmax)+b) + offset
    end
end

function Profiles(x::TabulatedTrento{A,B}, y::TabulatedTrento{A,B}, cent1::Integer, cent2::Integer; radius = range(0,30,100), norm_temp = 1, norm_coll = 1, xmax = 8, exp_tail = true) where {A,B}
    #temperature profile
    r, entropy_profile = get_profile(x, cent1, cent2; norm = norm_temp)
    
    temperature_profile = InverseFunction(x->pressure_derivative(x,Val(1),FluiduMEoS())).(entropy_profile) #.+0.0001 
    temperature_funct = LinearInterpolation(r, temperature_profile; extrapolation_bc=Flat())
    
    temp_exp = exponential_tail_pointlike.(Ref(temperature_funct), radius, Ref(xmax))
    temp_exp_funct = LinearInterpolation(radius, temp_exp; extrapolation_bc=Flat()) 
    
    #ncoll profile
    r, ncoll_profile = get_profile(y, cent1, cent2; norm = norm_coll) #.+0.0001 
    ncoll_funct = LinearInterpolation(r, ncoll_profile; extrapolation_bc=Flat()) 
    
    if exp_tail == true
        return temp_exp_funct, ncoll_funct
    else
        return temperature_funct, ncoll_funct
    end    
end


function Profiles_offset(x::TabulatedTrento{A,B}, y::TabulatedTrento{A,B}, cent1::Integer, cent2::Integer; radius = range(0,30,100), norm_temp = 1, norm_coll = 1, xmax = 8) where {A,B}
    #temperature profile
    r, entropy_profile = get_profile(x, cent1, cent2; norm = norm_temp)
    
    temperature_profile = InverseFunction(x->pressure_derivative(x,Val(1),FluiduMEoS())).(entropy_profile) #.+0.0001 
    temperature_funct = LinearInterpolation(r, temperature_profile; extrapolation_bc=Flat())
    
    
    #ncoll profile
    r, ncoll_profile = get_profile(y, cent1, cent2; norm = norm_coll) #.+0.0001 
    ncoll_funct = LinearInterpolation(r, ncoll_profile; extrapolation_bc=Flat()) 
    
    #return temp_exp_funct, ncoll_funct
    return temperature_funct, ncoll_funct    
end



#not used
function map_initial_profile(temperature, x::TabulatedTrento{A,B}, cent1::Integer, cent2::Integer; norm = 1, xmax = 8) where {A,B}
    r, profile = get_profile(x, cent1, cent2; norm = norm) 
    
    #interpolate the temperature which is itself a function of the entropy 
    temp_interpolated = LinearInterpolation(r, temperature; extrapolation_bc=Flat())
    #temp_interpolated = LinearInterpolation(r, funct_to_interpolate; extrapolation_bc=Flat()) 
    temp_exp = exponential_tail_pointlike.(Ref(temp_interpolated), r, Ref(xmax))  
    return LinearInterpolation(r, temp_exp; extrapolation_bc=Flat())  
end
#vectorializing now:
#exponential = exponential_tail_pointlike.((get_interpolate_profile(pallino2, 10, 20), ), vect, (xmax,))
#exponential = exponential_tail_pointlike.(Ref(get_interpolate_profile(pallino2, 10, 20)), vect, Ref(xmax))
