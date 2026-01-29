abstract type AbstractInitialCondition end

struct TabulatedData{A<:AbstractArray, B<:AbstractArray} <: AbstractInitialCondition
    radius::A
    profile::B
end

function TabulatedData(fn::A) where {A<:AbstractString}
    file = readdlm(fn,comment_char='#',comments=true)
    TabulatedData(file[:, 1], file[:, 2:end])
end   


function centrality_bin(x::TabulatedData{A,B}) where {A,B}
    n_columns = size(x.profile)[2]
    bin_size = div(100, n_columns)
end

"""
Get trento profiles when multiple columns corresponding to different centrality classes are given 
"""
function get_profile(x::TabulatedData{A,B}, cent1::Integer, cent2::Integer; norm = 1) where {A,B}
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

"""
Get trento profiles when only a single columns corresponding to a fixed centrality class is given 
"""
function get_profile(x::TabulatedData{A,B}; norm = 1) where {A,B}
    r = x.radius
    
    radius = r[sortperm(r)] #ensure that r is sorted ...

    mean_profile = vec(x.profile)
    mean_profile = mean_profile[sortperm(r)]  #... and mean_profile is sorted accordingly

    return radius, norm.*mean_profile
end

function exponential_tail_pointlike(function_profile,x;xmax = 8, offset = 0.015)
    a = ForwardDiff.derivative(function_profile, xmax) / (function_profile(xmax) - offset)
    b = log(function_profile(xmax) - offset) 
    if x-xmax < 0
        return function_profile(x)
    else  
        return exp(a*(x-xmax)+b) + offset
    end
end


function inverse_exponential_tail_pointlike(function_profile,x;xmax = 8, offset = 11)
    a = ForwardDiff.derivative(function_profile, xmax) / function_profile(xmax)
    b = log(-function_profile(xmax) + offset) 
    if x-xmax < 0
        return function_profile(x)
    else  
        return -exp(a*(x-xmax)+b) + offset
    end
end

function inverse_exponential_tail_pointlike_diff(function_profile, x; xmax, offset, h = 1e-6)
    # Use finite difference instead of ForwardDiff
    a = (function_profile(xmax + h) - function_profile(xmax - h)) / (2 * h) / (function_profile(xmax)-offset)
    b = log(-function_profile(xmax) + offset) 
    if x - xmax < 0
        return function_profile(x)
    else  
        return -exp(a*(x-xmax)+b) + offset
    end
end


function exponential_tail_pointlike(function_profile,x,y ;rmax = 8, offset = 0.015)
    theta = atan(y,x)
    xmax = rmax*cos(theta)
    ymax = rmax*sin(theta)
    a = ForwardDiff.gradient(v -> function_profile(v[1], v[2]), [xmax, ymax])/ function_profile(xmax,ymax)
    b = log(function_profile(rmax*cos(theta),rmax*sin(theta)) - offset) 
    if x^2+y^2 < rmax^2
        return function_profile(x,y)
    else  
        return exp(a[1]*(x-xmax)+a[2]*(y-ymax)+b) + offset
    end
end

"""
Initialize a profile when both the entropy and the density of binary collision profiles are given"""
function Profiles(x::TabulatedData{A,B}, y::TabulatedData{A,B}, cent1::Integer, cent2::Integer; radius, norm_x, norm_y, exp_tail = true,offset=0.005, xmax_temp = 8, xmax_ncoll = 5, offset_ncoll = 0.) where {A,B}
    #entropy profile
    r, entropy_profile = get_profile(x, cent1, cent2; norm = norm_x)
    
    temperature_profile = InverseFunction(x->pressure_derivative(x,Val(1),FluiduMEoS())).(entropy_profile)  
    temperature_funct = linear_interpolation(r, temperature_profile; extrapolation_bc=Flat())
    
    temp_exp = exponential_tail_pointlike.(Ref(temperature_funct), radius; xmax = xmax_temp, offset = offset)
   # temp_exp = exponential_tail_pointlike.(Ref(temperature_funct), radius; xmax = 8, offset = 0.01)
    temp_exp_funct = linear_interpolation(radius, temp_exp; extrapolation_bc=Flat()) 
    #ncoll profile
    r, ncoll_profile = get_profile(y, cent1, cent2; norm = norm_y)  
    ncoll_funct = linear_interpolation(r, ncoll_profile; extrapolation_bc=Flat()) 
   
    ncoll_exp = exponential_tail_pointlike.(Ref(ncoll_funct), radius;xmax = xmax_ncoll, offset= offset_ncoll)
        
    ncoll_exp_funct = linear_interpolation(radius, ncoll_exp; extrapolation_bc=Flat()) 
    
    if exp_tail == true
        return temp_exp_funct, ncoll_exp_funct
    else
        return temperature_funct, ncoll_funct
    end    
end


"""
Initialize when only one profile is given"""
function Profiles(x::TabulatedData{A,B}, cent1::Integer, cent2::Integer; radius, norm, xmax = 8, exp_tail = true, temperature_flag = true, offset = 0.01) where {A,B}
    #entopy profile
    r, entropy_profile = get_profile(x, cent1, cent2; norm = norm)
    
    if temperature_flag == true
        temperature_profile = InverseFunction(x->pressure_derivative(x,Val(1),FluiduMEoS())).(entropy_profile)  
        temperature_funct = linear_interpolation(r, temperature_profile; extrapolation_bc=Flat())
    else
        temperature_funct = linear_interpolation(r, entropy_profile; extrapolation_bc=Flat())
    end

    temp_exp = exponential_tail_pointlike.(Ref(temperature_funct), radius; xmax, offset)
    temp_exp_funct = linear_interpolation(radius, temp_exp; extrapolation_bc=Flat()) 
        
    if exp_tail == true
        return temp_exp_funct
    else
        return temperature_funct
    end    
end

"""
Initialize a profile when only the temperature profile is given"""
function Profiles(x::TabulatedData{A,B};radius = range(0,30,100), norm = 1, xmax = 8, exp_tail = true) where {A,B}
    #temperature profile
    r, temperature_profile = get_profile(x; norm = norm)
    temperature_funct = linear_interpolation(r, temperature_profile; extrapolation_bc=Flat())
        
    if exp_tail == true
    
        temp_exp = exponential_tail_pointlike.(Ref(temperature_funct), radius; xmax)
        temp_exp_funct = linear_interpolation(radius, temp_exp; extrapolation_bc=Flat()) 
        return temp_exp_funct
    else
        return temperature_funct
    end   
end


