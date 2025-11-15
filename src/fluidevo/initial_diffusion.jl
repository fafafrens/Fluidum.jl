"""
Calculate the density of binary collisions after free streaming is applied"""
function ncoll_fs(ncoll_profile,x,y,px,py; initial_params, m = 1.5) 
    tau0 = initial_params.tau0
    ptau = sqrt(px^2 + py^2 + m^2) 
    x1 = x - px*tau0/ptau
    y1 = y - py*tau0/ptau
    return ncoll_profile(sqrt(x1^2+y1^2))
end 

"""
Take fonll production cross section"""
function dσ_fonll(fonll_profile,px,py)  
    return fonll_profile((sqrt(px^2+py^2)))
end

"""
Consider the numerical integration of differential cross-section multiplied by ncoll"""
function A_prime(fonll_profile,ncoll_profile,x,y;initial_params, m = 1.5) 
    fonll_function(px,py) = dσ_fonll(fonll_profile,px,py) 
    ncoll_shift(px,py) = ncoll_fs(ncoll_profile, x,y,px,py;initial_params)
    hcubature(p->(fonll_function(p[1],p[2])*ncoll_shift(p[1],p[2])), rtol=0.001, [-20/sqrt(2), -20/sqrt(2)], [20/sqrt(2), 20/sqrt(2)])[1]
end


"""
Calculate equilibrium production cross section normalized to reproduce A"""
function dσ_eq_new(T,fonll_profile, ncoll_profile,x,y,px,py; initial_params, m=1.5)
    mt = sqrt(px^2 + py^2 + m^2)
    A = A_prime(fonll_profile,ncoll_profile,x,y;initial_params) 
    ncoll_shift = ncoll_fs(ncoll_profile, x,y,px,py;initial_params)      
    dσ_equil = 1/(2*pi)*(mt*besselk(1,mt/T))/(m^2*T*besselk(2,m/T))*A/ncoll_shift #equilibrium cross section
    return dσ_equil
end


"""
Calculate the charm density using the equilibrium distribution function, after free streaming is applied"""
function density_fs_equilibrium(T,fonll_profile,ncoll_profile,r; initial_params, m = 1.5) 
    tau0 = initial_params.tau0

    eq_interp(x,y,px,py) = dσ_eq_new(T,fonll_profile, ncoll_profile,x,y,px,py; m,initial_params)
    ncoll_shift(x,y,px,py) = ncoll_fs(ncoll_profile,x,y,px,py;initial_params) 
    density(x,y) = 1/(tau0)*hcubature(p->(ncoll_shift(x,y,p[1],p[2])*eq_interp(x,y,p[1],p[2])), rtol=0.0001, [-20/sqrt(2), -20/sqrt(2)], [20/sqrt(2), 20/sqrt(2)])[1]
    density_polar = density(r,0)
    return density_polar
end 

"""
Calculate the charm density using the total distribution function, after free streaming is applied (should match with density_fs_equilibrium)"""
function density_fs_total(fonll_profile,ncoll_profile,r; initial_params) 
    tau0 = initial_params.tau0
    
    fonll_interp(px,py) = dσ_fonll(fonll_profile,px,py)
    ncoll_shift(x,y,px,py) = ncoll_fs(ncoll_profile,x,y,px,py;initial_params) 
    density(x,y) = 1/(tau0)*hcubature(p->(ncoll_shift(x,y,p[1],p[2])*fonll_interp(p[1],p[2])), rtol=0.001, [-20/sqrt(2), -20/sqrt(2)], [20/sqrt(2), 20/sqrt(2)])[1]
    density_polar = density(r,0)
    return density_polar
end 

"""
Calculate the diffusion current knowing that only total distribution function contributes"""
function nur(T, fonll_profile,ncoll_profile,r; initial_params, m = 1.5)
    tau0 = initial_params.tau0

    ncoll_shift(x,y,px,py) = ncoll_fs(ncoll_profile,x,y,px,py;initial_params) 
    fonll_interp(px,py) = dσ_fonll(fonll_profile,px,py)
    eq_interp(x,y,px,py) = dσ_eq_new(T, fonll_profile, ncoll_profile,x,y,px,py; initial_params, m)

    fonll_integral(x,y) = hcubature(p->(p[1]/sqrt(m^2+p[1]^2+p[2]^2)*ncoll_shift(x,y,p[1],p[2])*fonll_interp(p[1],p[2])), rtol=0.001, [-20/sqrt(2), -20/sqrt(2)], [20/sqrt(2), 20/sqrt(2)])[1]
    #eq_integral(x,y) = hcubature(p->(p[1]/sqrt(m^2+p[1]^2+p[2]^2)*ncoll_shift(x,y,p[1],p[2])*eq_interp(x,y,p[1],p[2])), rtol=0.001, [-20/sqrt(2), -20/sqrt(2)], [20/sqrt(2), 20/sqrt(2)])[1]
    #nux(x,y) = 1/(tau0)*(fonll_integral(x,y) - eq_integral(x,y))
    nux(x,y) = 1/(tau0)*fonll_integral(x,y)
    nur_polar = nux(r,0)
    return nur_polar
end

"""
Calculate the equilibrium contribution to the diffusion current (should be zero)"""
function nur_eq(T, fonll_profile,ncoll_profile,r; initial_params, m = 1.5)
    tau0 = initial_params.tau0

    ncoll_shift(x,y,px,py) = ncoll_fs(ncoll_profile,x,y,px,py;initial_params) 
    eq_interp(x,y,px,py) = dσ_eq_new(T, fonll_profile, ncoll_profile,x,y,px,py; initial_params, m)

    eq_integral(x,y) = hcubature(p->(p[1]/sqrt(m^2+p[1]^2+p[2]^2)*ncoll_shift(x,y,p[1],p[2])*eq_interp(x,y,p[1],p[2])), rtol=0.01, [-20/sqrt(2), -20/sqrt(2)], [20/sqrt(2), 20/sqrt(2)])
    nux_eq(x,y) = 1/(tau0).*eq_integral(x,y)
    nur_polar = nux_eq(r,0)

    return nur_polar
end


