#using SciMLBase #for ODEProblem
#using OrdinaryDiffEq #for ODEProblem

"""
    integral_cauchy(resultNofo,time,grid)

    resultNofo is the result of oneshoot
    grid = discretization.grid
    
    Returns values of the conserved charge at t = time

"""
particle_list=pwd()*"/../src/kernels/particles_beauty_pythia8_list.txt"

function integral_cauchy(resultNofo,time,grid,eos)
    mu(time)=LinearInterpolation([grid[i][1] for i in 2:lastindex(grid)-1],resultNofo(time)[6,2:lastindex(grid)-1]; extrapolation_bc=Flat())
    t(time)=LinearInterpolation([grid[i][1] for i in 2:lastindex(grid)-1],resultNofo(time)[1,2:lastindex(grid)-1]; extrapolation_bc=Flat())
    nu(time)=LinearInterpolation([grid[i][1] for i in 2:lastindex(grid)-1],resultNofo(time)[7,2:lastindex(grid)-1]; extrapolation_bc=Flat())
    u(time)=LinearInterpolation([grid[i][1] for i in 2:lastindex(grid)-1],resultNofo(time)[2,2:lastindex(grid)-1]; extrapolation_bc=Flat())
    ut(time,x) = sqrt(1+u(time)(x)*u(time)(x)) 
    nt(time,x) = u(time)(x)*nu(time)(x)/ut(time,x)
    density(time,x) = thermodynamic(t(time)(x),mu(time)(x),eos.hadron_list; ccbar = 1.4).pressure
    #@show quadgk(x->density(time,x),0,30)
    #density(time,x)=free_charm(t(time)(x),mu(time)(x),Heavy_Quark(mass=4.75,hadron_list=HadronResonaceGas(name_file=particle_list,Maxmass=11.,Minmass=4.0)))[1]
    
    return quadgk(x->2*pi*x*time*(density(time,x)* ut(time,x) +nt(time,x)),grid[2][1],grid[end-1][1],rtol=0.00001)

end

"""
    test_integral_cauchy(resultNofo,grid,tspan;dt=0.5)

    resultNofo is the result of oneshoot
    grid = discretization.grid

    Returns values of the conserved charge at each time in tspan; timestep can be specified with dt.

"""
function test_integral_cauchy(resultNofo,grid,tspan,eos;dt=0.5)
    current = []
    for time in tspan[1]:dt:tspan[2]
        push!(current,integral_cauchy(resultNofo,time,grid,eos))
    end
    return current
end
"""
    fo_integrand(alpha,x,phi)
    
    alpha is the freeze-out parameter
    x are the t,r coordinates as function of alpha
    phi are the fields as function of alpha
    
    Computes the conserved current on the freeze out Surface

"""

function fo_integrand(alpha,x,phi,eos)
    t,r= x(alpha)
    dta,dra=jacobian(x,alpha)
    T,ur,pi_phi,pi_eta,pi_b,μ,ν=phi(alpha)
    ur = abs(ur)
    #n = federica(T,μ,Heavy_Quark())[1]
    thermo = thermodynamic(T,μ,eos.hadron_list)
    n=thermo.pressure
    ut = sqrt(1+ur*ur) 
    nt = ur*ν/ut
    factor = t*r*2*pi
    #nt=0
    #nr=0
    return (-(n*ut+nt)*dra+ (n*ur+ν)*dta)*factor
    #return (n*ut)*dra
end

"""
    fo_integral(fo::FreezeOutResult{A,B};rtol=0.00001)

    Computes integral of the conserved current on the freeze-out surface. It should be compatible with the integral_cauchy at the latest freezeout time.
"""

function fo_integral(fo::FreezeOutResult{A,B},eos;rtol=0.00001) where {A<:SplineInterp,B<:SplineInterp}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)
    return quadgk(alpha->fo_integrand(alpha,x,phi,eos),lb...,rb...,rtol=rtol)
end
"""
function freeze_out_routine(oned_visc_hydro_discrete,matrxi1d_visc_HQ!,fluidpropery,phi,tspan)
returns freeze-out coordinates and fields values at freeze out
"""
function freeze_out_routine(oned_visc_hydro_discrete,matrxi1d_visc_HQ!,fluidpropery,phi,tspan;Tfo=0.1565)
    result=isosurface(oned_visc_hydro_discrete,matrxi1d_visc_HQ!,fluidpropery,phi,tspan,:temperature,Tfo);
    cha=Chart(Surface(result[:surface]),(t,x)->SVector{1}(atan(t,x)))
    return freezeout_interpolation(cha,sort_index=1)
end



