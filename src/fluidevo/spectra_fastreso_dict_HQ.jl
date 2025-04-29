
@inline function cilindrical_thermal_spectra(q,r,t,dra,dta,T,ν,μ,norm,K1eq,K2eq,K1diff,K2diff; delta_f = false)
    

    fmGeV = 5.068
    #mt=sqrt(m^2+pt^2)
    factor=1/(2*pi)^3*t*r
    
    r_factor=-factor* dra
    t_factor=factor* dta
    
    result = (r_factor.*(K1eq.+K1diff*ν*q./norm)+t_factor*(K2eq.+K2diff*ν*q./norm)).*exp(μ)
    
    if delta_f == false
        ν = 0
        n = 1
       result = (r_factor*(K1eq+K1diff*ν/(T*n))+t_factor*(K2eq+K2diff*ν/(T*n)))*exp(μ)    
    end
       
    return result*fmGeV^3
end


@inline function _pointwise_spectra(pt,alpha,x::A,phi::B,part::particle_attribute{S,R,U,V};decays, ccbar = 2.76, delta_f = false) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    
    fact = besseli(1, ccbar/2)/besseli(0, ccbar/2)

    t,r= x(alpha)
    dta,dra=jacobian(x,alpha)
    T,ur,pi_phi,pi_eta,pi_b,μ,ν=phi(alpha) #prendono gli argomenti da phi, che li prende da fo che viene calcolato altrove
    
    eos = Heavy_Quark()
    norm = normalization(T,μ,eos)
    
    
    if (ur > 5) 
        @warn string("Radial velocity out of fastreso limits")       
    end 
    
    q = 1
    if part.name == "Dc2007zer" || part.name == "Dc2010plu" #D* ha degenerazione 3 perche ha total angular momentum = 3. Per tutti gli altri ho spin 0 e non va preso in considerazione 
        deg = 3
    else deg = 1
    end

    if part.name == "jp3096zer"
        fact = 1
        q = 2
    end 
    μ = q*μ 

    if(decays==true)
        kernel = part.total_kernel_ext
    else 
        kernel = part.thermal_kernel_ext
    end
    K1eq=kernel.K1eq(T,pt,ur)
    K2eq=kernel.K2eq(T,pt,ur)
    K1diff=kernel.K1diff(T,pt,ur)
    K2diff=kernel.K2diff(T,pt,ur)
    cilindrical_thermal_spectra(q,r,t,dra,dta,T,ν,μ,norm,K1eq,K2eq,K1diff,K2diff;delta_f=delta_f)*deg*fact
end

#spectra in a single pt point
function spectra(pt::C,fo::FreezeOutResult{A,B},part::particle_attribute{S,R,U,V};rtol=1000*sqrt(eps()),decays=true,ccbar = 2.76) where {C<:Number,A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    x,phi=fo   
    lb=leftbounds(x)
    rb=rightbounds(x)

    #numerical integration of function of alpha, from the lb of alpha to the rb of alpha (over all the fo surface)
    quadgk(alpha->_pointwise_spectra(pt,alpha,x,phi,part;decays,ccbar,delta_f),lb...,rb...,rtol=rtol)

end

#spectra points between max, min with given step (even spacing)
function spectra(fo::FreezeOutResult{A,B},part::particle_attribute{S,R,U,V};pt_min=0.,pt_max=10.0,step=100,rtol=1000*sqrt(eps()),decays=true,rightbound=(100,),ccbar = 2.76,delta_f=false) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    x,phi=fo
    
    #T = [fo.fields.a[i][1] for i in 1:100]
    #μ = [fo.fields.a[i][6] for i in 1:100]

    #eos = Heavy_Quark()
    #norm = normalization.(T,μ,Ref(eos))
    #@show T
    #@show typeof(T)
    #norm = 1
    
    lb=leftbounds(x)
    rb=rightbounds(x)
    rb=min.(rb,rightbound)
    buff=alloc_segbuf(Float64, eltype(lb),Float64 ;size=1)
    
    if delta_f == true print("delta_f applied \n")  
    else print("delta_f not applied \n")
    end 

    [quadgk(alpha->_pointwise_spectra(pt,alpha,x,phi,part;decays,ccbar,delta_f),lb...,rb...;segbuf=buff,rtol=rtol) for pt in range(pt_min,pt_max,step) ] 

end

#spectra points in a given range (can also be uneven spacing)
function spectra(fo::FreezeOutResult{A,B},part::particle_attribute{S,R,U,V}, pt_range::Vector{Float64};rtol=1000*sqrt(eps()),decays=true,rightbound=(100,),ccbar = 2.76) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)
    rb=min.(rb,rightbound)
    buff=alloc_segbuf(Float64, eltype(lb),Float64 ;size=1)
    
    [quadgk(alpha->_pointwise_spectra(pt,alpha,x,phi,part;decays,ccbar),lb...,rb...;segbuf=buff,rtol=rtol) for pt in pt_range ] 

end


function spectra(pt::C,fo::FreezeOutResult{A,B},part::particle_attribute{S,R,U,V};rtol=1000*sqrt(eps()),decays=true) where {C<:AbstractVector,A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)
    buff=alloc_segbuf(Float64, eltype(lb),Float64;size=1)

    [quadgk(alpha->_pointwise_spectra(i,alpha,x,phi,part;decays),lb...,rb...,segbuf=buff,rtol=rtol) for i in pt ] 

end

#spectra points over pt range
function spectra(m::Number,fo::FreezeOutResult{A,B};pt_min=0.,pt_max=8.0,step=100) where {A<:SplineInterp,B<:SplineInterp}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)
    buff=alloc_segbuf(Float64, eltype(lb),Float64 ;size=1)
    
    [quadgk(alpha->_pointwise_spectra(pt,m,alpha,x,phi),lb...,rb...;segbuf=buff) for pt in range(pt_min,pt_max,step) ] 

end
function multiplicity(pt,spectra)
    f = [spectra[i][1]*pt[i] for i in 1:length(pt)]
    2.0*π*NumericalIntegration.integrate(pt,f)
end


function multiplicity(fo::FreezeOutResult{A,B},part::particle_attribute{S,R,U,V};rtol=1000*sqrt(eps()),decays=true, delta_f = false, rightbound=100,pt_min=0.,pt_max=10.0) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)
    rb = min.(rb,(rightbound,))
    hcubature( b ->2.0*π *_pointwise_spectra(b[1],b[2],x,phi,part;decays,delta_f)*b[1],(pt_min,lb...),(pt_max,rb...);rtol=rtol)
    
end

function _pointwise_spectra_internal(pt,m,alpha,x::A,phi::B;deg=1) where {A<:SplineInterp,B<:SplineInterp}
    t,r= x(alpha)
    dta,dra=jacobian(x,alpha)
    T,ur,pi_phi,pi_eta,pi_b,μ,ν=phi(alpha)
    internal_thermal_spectra(pt,m,r,t,dra,dta,ur,T,μ,ν,deg=deg)
end
function internal_thermal_spectra(pt,m,r,t,dra,dta,ur,T,μ,ν;deg=1)
    fmGeV = 5.068
    mt=sqrt(m^2+pt^2)
    factor=1/(2*pi^2)*t*r
    karg=mt*sqrt(1+ur^2)/T
    iarg=pt*ur/T
    #n = federica(T,μ,Heavy_Quark())[1]
    #n = 1
    ut = sqrt(1+ur*ur) 
    
    r_factor=-factor* dra*mt
    t_factor=factor* dta*pt

    i0=besseli0(iarg)
    i1=besseli1(iarg)
    i2=besseli(2,iarg)
    i2min=besseli(-2,iarg)
    k0=besselk0(karg)
    k1=besselk1(karg)
    k2=besselk(2,karg)
    k2min=besselk(-2,karg)
    rcc = 1/4*(k2min+2*k0+k2)
    acc = 1/4*(i2min+2*i0+i2)
    #ν=-ν   #the correction has a minus sign in front of the 1/P_hq
    ν = 0 #commented on 23.10
    #n = thermodynamic(T,μ,eos.hadron_list).pressure
    n=1

    result= (r_factor*(k1*i0-i1*k1*pt*ν/n/T+mt*ν/n/T*ur/ut*i0*rcc)+
    t_factor*(k0*i1-pt*k0*ν/n/T*acc+mt*ν/n/T*ur/ut*i1*k1))*exp(μ)
    #result= (r_factor*(k1*i0)+t_factor*(k0*i1))*exp(μ)
    #corr = (r_factor*(-i1*k1*pt*ν/n/T+mt*ν/n/T*ur/ut*i0*rcc)+
    #t_factor*(-pt*k0*ν/n/T*acc+mt*ν/n/T*ur/ut*i1*k1))*exp(μ)
    #corr = (r_factor*(-i1*k1*pt*ν/n/T+mt*ν/n/T*ur/ut*i0*rcc))*exp(μ)
    #corr = (t_factor*(-pt*k0*ν/n/T*acc+mt*ν/n/T*ur/ut*i1*k1))*exp(μ)
    
    return result*fmGeV^3*deg
end 

function spectra_internal(m::Number,fo::FreezeOutResult{A,B};pt_min=0.,pt_max=8.0,step=100,deg=1) where {A<:SplineInterp,B<:SplineInterp}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)
    buff=alloc_segbuf(Float64, eltype(lb),Float64 ;size=1)
    
    [quadgk(alpha->_pointwise_spectra_internal(pt,m,alpha,x,phi,deg=deg),lb...,rb...;segbuf=buff) for pt in range(pt_min,pt_max,step) ] 

end



function multiplicity(m::Number,fo::FreezeOutResult{A,B};rtol=sqrt(eps()),pt_min=0,pt_max = 10) where {A<:SplineInterp,B<:SplineInterp}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)

    #hcubature refers to a high-dimensional integration routine provided by the HCubature.jl package. It is used for numerical integration over hypercubes in multiple dimensions.
    hcubature( b ->2.0*π *_pointwise_spectra_internal(b[1],m,b[2],x,phi)*b[1],(pt_min,lb...),(pt_max,rb...);rtol=rtol)
    
end
#WITH DISSIPATIVE CORRECTIONS
struct prefactors{T}
    piphiphi::T
    pietaeta::T
    bulkfactor::T
end

function prefactors(T,pi_phi,pi_eta,pi_b,fluidpropery)
    #ur=field[1], piϕϕ->field[2], piηη->field[3], πbulk -> field[4], νr -> field[5], μ=field[6]
    dtp=pressure_derivative(T,Val(1),fluidpropery)
    #dtp=1.68
    dtdtp=pressure_derivative(T,Val{2}(),fluidpropery.eos)
    
    
    #non-zero bulk case
    tauB=τ_bulk(T,dtp,dtdtp,fluidpropery.bulk)
    zeta=bulk_viscosity(T,dtp,fluidpropery.bulk)
    

    #eta=viscosity(temperature,dtp,fluidpropery.shear)
    
    bulkfactor = pi_b*tauB/zeta
    if (zeta == 0.0)
        bulkfactor=0.0
    end
    pietaeta = pi_eta/(2*(T*dtp)*T^2)
    piphiphi = pi_phi/(2*(T*dtp)*T^2)
    
    return prefactors(piphiphi, pietaeta, bulkfactor)
end

@inline function cilindrical_thermal_spectra(pt,m,r,t,dra,dta,ur,T,pi_phi,pi_eta,pi_b,K1eq,K2eq,K1piphi,K2piphi,K1pieta,K2pieta,K1bulk,K2bulk,fluidpropery)

    @unpack piphiphi,pietaeta,bulkfactor = prefactors(T,pi_phi,pi_eta,pi_b,fluidpropery)
    fmGeV = 5.068
    #mt=sqrt(m^2+pt^2)
    factor=1/(2*pi)^3*t*r
    
    r_factor=-factor* dra
    t_factor=factor* dta
    
    result= r_factor*(K1eq+pietaeta*K1pieta+piphiphi*K1piphi+
    bulkfactor*K1bulk)+t_factor*(K2eq+pietaeta*K2pieta+piphiphi*K2piphi+bulkfactor*K2bulk) 
    #print("all")
    #result= r_factor*(K1eq)+t_factor*(K2eq)
    
    return result*fmGeV^3
end

#used for lf 
@inline function _pointwise_spectra(pt,alpha,x::A,phi::B,part::particle_attribute{S,R,U,V},fluidpropery::FluidProperties{C,D,E,F};decays) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V,C,D,E,F}
    t,r= x(alpha)
    dta,dra=jacobian(x,alpha)
    T,ur,pi_phi,pi_eta,pi_b=phi(alpha)
    m = part.mass
    #@show ur 
    if (ur > 3) 
        @show ur 
        @warn string("Radial velocity out of fastreso limits")       
    end

    if(decays==true)
        kernel = part.total_kernel_ext
    else 
        kernel = part.thermal_kernel_ext
    end
    K1eq=kernel.K1eq(T,pt,ur)
    K2eq=kernel.K2eq(T,pt,ur)
    K1piphi=kernel.K1piphi(T,pt,ur)
    K2piphi=kernel.K2piphi(T,pt,ur)
    K1pieta=kernel.K1pieta(T,pt,ur)
    K2pieta=kernel.K2pieta(T,pt,ur)
    K1bulk=kernel.K1bulk(T,pt,ur)
    K2bulk=kernel.K2bulk(T,pt,ur)
    #print("fixed temperature")
    # K1eq=kernel.K1eq(pt,ur)
    # K2eq=kernel.K2eq(pt,ur)
    # K1piphi=kernel.K1piphi(pt,ur)
    # K2piphi=kernel.K2piphi(pt,ur)
    # K1pieta=kernel.K1pieta(pt,ur)
    # K2pieta=kernel.K2pieta(pt,ur)
    # K1bulk=kernel.K1bulk(pt,ur)
    # K2bulk=kernel.K2bulk(pt,ur)
    
    cilindrical_thermal_spectra(pt,m,r,t,dra,dta,ur,T,pi_phi,pi_eta,pi_b,K1eq,K2eq,K1piphi,K2piphi,K1pieta,K2pieta,K1bulk,K2bulk,fluidpropery)
end



#spectra in a single pt point
function spectra(pt::C,fo::FreezeOutResult{A,B},part::particle_attribute{S,R,U,V},fluidpropery;rtol=1000*sqrt(eps()),decays=true) where {C<:Number,A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)

    #numerical integration of function of alpha, from the lb of alpha to the rb of alpha (over all the fo surface)
    #returns a pair (integral, error)
    quadgk(alpha->_pointwise_spectra(pt,alpha,x,phi,part,fluidpropery;decays),lb...,rb...,rtol=rtol)

end






#spectra points over pt range
function spectra(fo::FreezeOutResult{A,B},part::particle_attribute{S,R,U,V},fluidpropery;pt_min=0.,pt_max=10.0,step=100,rtol=1000*sqrt(eps()),decays=true) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)
    buff=alloc_segbuf(Float64, eltype(lb),Float64 ;size=1)
    
    [quadgk(alpha->_pointwise_spectra(pt,alpha,x,phi,part,fluidpropery;decays),lb...,rb...;segbuf=buff,rtol=rtol) for pt in range(pt_min,pt_max,step) ] 

end

#spectra points over pt range with uneven spacing
function spectra_lf(fo::FreezeOutResult{A,B},part::particle_attribute{S,R,U,V},fluidpropery, pt_range::Vector{Float64};rtol=1000*sqrt(eps()),decays=true) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)
    buff=alloc_segbuf(Float64, eltype(lb),Float64 ;size=1)
    
    [quadgk(alpha->_pointwise_spectra(pt,alpha,x,phi,part,fluidpropery;decays),lb...,rb...;segbuf=buff,rtol=rtol) for pt in pt_range ] 

end

function spectra(pt::C,fo::FreezeOutResult{A,B},part::particle_attribute{S,R,U,V},fluidpropery;rtol=1000*sqrt(eps()),decays=true) where {C<:AbstractVector,A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)
    buff=alloc_segbuf(Float64, eltype(lb),Float64;size=1)

    [quadgk(alpha->_pointwise_spectra(i,alpha,x,phi,part,fluidpropery;decays),lb...,rb...,segbuf=buff,rtol=rtol) for i in pt ] 

end

function multiplicity_lf(fo::FreezeOutResult{A,B},part::particle_attribute{S,R,U,V},fluidpropery;rtol=1000*sqrt(eps()),decays=true, rightbound=100,pt_min=0.,pt_max=10.0) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)
    rb = min.(rb,(rightbound,))
    deg = part.degeneracy
    hcubature( b ->2.0*π *_pointwise_spectra(b[1],b[2],x,phi,part,fluidpropery;decays)*b[1],(pt_min,lb...),(pt_max,rb...);rtol=rtol)
    
end

