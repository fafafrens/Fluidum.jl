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
    a = ForwardDiff.derivative(function_profile, xmax) / function_profile(xmax)
    b = log(function_profile(xmax) - offset) 
    if x-xmax < 0
        return function_profile(x)
    else  
        return exp(a*(x-xmax)+b) + offset
    end
end


"""
Initialize a profile when both the entropy and the density of binary collision profiles are given"""
function Profiles(x::TabulatedData{A,B}, y::TabulatedData{A,B}, cent1::Integer, cent2::Integer; radius = range(0,30,100), norm_temp = 1, norm_coll = 1, exp_tail = true, offset = 0) where {A,B}
    #entropy profile
    r, entropy_profile = get_profile(x, cent1, cent2; norm = norm_temp)
    
    temperature_profile = InverseFunction(x->pressure_derivative(x,Val(1),FluiduMEoS())).(entropy_profile)  
    temperature_funct = linear_interpolation(r, temperature_profile; extrapolation_bc=Flat())
    
    temp_exp = exponential_tail_pointlike.(Ref(temperature_funct), radius; xmax = 8, offset = 0.005)
   # temp_exp = exponential_tail_pointlike.(Ref(temperature_funct), radius; xmax = 8, offset = 0.01)
    
    temp_exp_funct = linear_interpolation(radius, temp_exp; extrapolation_bc=Flat()) 
    
    #ncoll profile
    r, ncoll_profile = get_profile(y, cent1, cent2; norm = norm_coll)  
    ncoll_funct = linear_interpolation(r, ncoll_profile; extrapolation_bc=Flat()) 
    
    ncoll_exp = exponential_tail_pointlike.(Ref(ncoll_funct), radius;xmax = 5, offset)
        
    ncoll_exp_funct = linear_interpolation(radius, ncoll_exp; extrapolation_bc=Flat()) 
    
    if exp_tail == true
        return temp_exp_funct, ncoll_exp_funct
    else
        return temperature_funct, ncoll_funct
    end    
end


"""
Initialize when only one profile is given"""
function Profiles(x::TabulatedData{A,B}, cent1::Integer, cent2::Integer; radius = range(0,30,100), norm = 1, xmax = 8, exp_tail = true, temperature_flag = false, offset = 0) where {A,B}
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
function Profiles(x::TabulatedData{A,B};radius = range(0,30,100), norm_ = 1, xmax = 8, exp_tail = true) where {A,B}
    #temperature profile
    r, temperature_profile = get_profile(x; norm = norm_)
    temperature_funct = linear_interpolation(r, temperature_profile; extrapolation_bc=Flat())
        
    if exp_tail == true
    
        temp_exp = exponential_tail_pointlike.(Ref(temperature_funct), radius; xmax)
        temp_exp_funct = linear_interpolation(radius, temp_exp; extrapolation_bc=Flat()) 
        return temp_exp_funct
    else
        return temperature_funct
    end   
end

"""
Take fonll function"""
function fonll(fonll_profile,px,py; m = 1.5)  
    return fonll_profile((sqrt(px^2+py^2)))/(sqrt(px^2+py^2+m^2))
end

function dσ_eq(T,pt; dσ_QQdy, m=1.5)         
    dσ_equil(pt) = exp(-sqrt(pt^2+m^2)/T)/(2*pi*T*(m^2+2*m*T+2*T^2)*exp(-m/T))*dσ_QQdy #equilibrium cross section
    return dσ_equil(pt)
end


"""
Calculate the density of binary collisions after free streaming is applied"""
function ncoll_fs(ncoll_profile,x,y,px,py; tau0 = 0.4, m = 1.5) 
    ptau = sqrt(px^2 + py^2 + m^2) 
    x1 = x - px*tau0/ptau
    y1 = y - py*tau0/ptau
    return ncoll_profile(sqrt(x1^2+y1^2))
end 


"""
Calculate the charm density after free streaming is applied"""

function density_fs(T,ncoll_profile,x,y;dσ_QQdy, tau0, m = 1.5) 
    eq_interp(px,py) = dσ_eq(T, sqrt(px^2+py^2);dσ_QQdy, m)
    ncoll_shift(px,py) = ncoll_fs(ncoll_profile,x,y,px,py;tau0) 
    density = 1/(dσ_QQdy)*hcubature(p->(sqrt(m^2+p[1]^2+p[2]^2)*ncoll_shift(p[1],p[2])*eq_interp(p[1],p[2])), rtol=0.0000001, [-20/sqrt(2), -20/sqrt(2)], [20/sqrt(2), 20/sqrt(2)])[1]
    return density
end 


function density_polar(T,ncoll_profile,r;dσ_QQdy, tau0 = 0.4, m = 1.5) 
    density_fs_(x,y) = density_fs(T,ncoll_profile,x,y;dσ_QQdy, tau0, m) 
    density_polar = density_fs_(r,0) 
    return density_polar
end

function nux(T,ncoll_profile, fonll_profile,x,y; tau0, dσ_QQdy, m = 1.5)
    ncoll_shift(px,py) = ncoll_fs(ncoll_profile,x,y,px,py;tau0) 
    fonll_interp(px,py) = fonll(fonll_profile,px,py; m)
    eq_interp(px,py) = dσ_eq(T, sqrt(px^2+py^2);dσ_QQdy)

    fonll_integral = hcubature(p->(p[1]*ncoll_shift(p[1],p[2])*fonll_interp(p[1],p[2])), rtol=0.001, [-20/sqrt(2), -20/sqrt(2)], [20/sqrt(2), 20/sqrt(2)])[1]
    eq_integral = hcubature(p->(p[1]*ncoll_shift(p[1],p[2])*eq_interp(p[1],p[2])), rtol=0.001, [-20/sqrt(2), -20/sqrt(2)], [20/sqrt(2), 20/sqrt(2)])[1]
    nux = 1/(dσ_QQdy)*(fonll_integral - eq_integral)
    return nux
end


"""
Quality check: is the equilibrium cross section normalized to 1?"""
function dσ_eq_norm(T,pt; dσ_QQdy, m=1.5)         
    dσ_equil(pt) = exp(-sqrt(pt^2+m^2)/T)/(2*pi*T*(m^2+2*m*T+2*T^2)*exp(-m/T))*dσ_QQdy #equilibrium cross section
    normalization =  quadgk(x->2*pi*x*sqrt(x^2+m^2)*exp(-sqrt(x^2+m^2)/T),0.1,20,rtol=0.00001)[1]/(2*pi*T*(m^2+2*m*T+2*T^2)*exp(-m/T))  
    
    @show normalization 
    
    dσ_QQdy_calculated = quadgk(x->2*pi*x*sqrt(x^2+m^2)*dσ_equil(x),0.1,20,rtol=0.00001)[1]
    @show dσ_QQdy_calculated
    return dσ_equil(pt)
end




"""
Quality check: interpolation of the ncoll and of fonll"""
function interpp(data::TabulatedData{A,B}; cent1=0, cent2=10) where {A,B}
    r, ncoll_trento = get_profile(data, cent1, cent2; norm = 1) #ncoll profile
    ncoll_funct = linear_interpolation(r, ncoll_trento; extrapolation_bc=Flat())
    return ncoll_funct
end 

function dσ_fonll_interp(x::TabulatedData{A,B}) where {A,B}         
    dσ_fonll_interpol = LinearInterpolation(x.radius, vec(x.profile)*1e-10; extrapolation_bc=Flat()) #interpolated FONLL cross section
    return dσ_fonll_interpol
end

