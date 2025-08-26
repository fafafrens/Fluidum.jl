function cilindrical_thermal_spectra(q,pt,m,r,t,dra,dta,ur,T,ν,μ,eos; delta_f = false)
    fmGeV = 5.068
    mt=sqrt(m^2+pt^2)
    nut = ν*ur/sqrt(1+ur^2)
    factor=1/(2*pi^2)*t*r
    karg=mt*sqrt(1+ur^2)/T
    iarg=pt*ur/T

    r_factor=-dra*mt
    t_factor=dta*pt

    i0=besseli0(iarg)
    i1=besseli1(iarg)
    k0=besselk0(karg)
    k1=besselk1(karg)
    eq = (r_factor*k1*i0+t_factor*k0*i1)*factor
    result = eq*exp(μ)

    if delta_f == true
        factor1=dra*mt^2*nut
        factor2=-dra*mt*pt*ν
        factor3=-dta*mt*pt*nut
        factor4=dta*pt^2*ν

        idiff = 1/2*(besseli0(iarg) + besseli(2,iarg))
        kdiff = 1/2*(besselk0(karg) + besselk(2,karg))
        
        norm = normalization(T,μ,eos)
        diff = (factor1*kdiff*i0+factor2*k1*i1+factor3*k1*i1+factor4*k0*idiff)*factor*(q/(T*norm))
        result = (eq+diff)*exp(μ)
    end

    #result= (r_factor*k1*i0+t_factor*k0*i1)*exp(μ)
    return result*fmGeV^3
end 


function cilindrical_thermal_spectra_exponential(q,pt,m,r,t,dra,dta,ur,T,ν,μ,eos; delta_f = false)
    fmGeV = 5.068
    mt=sqrt(m^2+pt^2)
    nut = ν*ur/sqrt(1+ur^2)
    factor=1/(2*pi^2)*t*r
    norm = normalization(T,μ,eos)
    karg=mt/T*(sqrt(1+ur^2) + nut*q/norm)
    iarg=pt/T*(ur + ν*q/norm)

    r_factor=-dra*mt
    t_factor=dta*pt
    
    i0=besseli0(iarg)
    i1=besseli1(iarg)
        
    k0=besselk0(karg)
    k1=besselk1(karg)
    result= (r_factor*k1*i0+t_factor*k0*i1)*factor*exp(μ)
    
    return result*fmGeV^3
end 


function _pointwise_spectra_analytic(pt,alpha,x::A,phi::B,part::particle_attribute{S,R,U,V},eos;delta_f) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    t,r= x(alpha)
    dta,dra=jacobian(x,alpha)
    T,ur,pi_phi,pi_eta,pi_b,μ,ν=phi(alpha)
    m = part.mass
    
    fact = besseli(1, eos.hadron_list.ccbar/2)/besseli(0, eos.hadron_list.ccbar/2)
    q = 1
    if part.name == "Dc2007zer" || part.name == "Dc2010plu"  
        deg = 3
    else deg = 1
    end

    if part.name == "jp3096zer"
        fact = 1
        q = 2
    end 
    μ = q*μ 

    cilindrical_thermal_spectra_exponential(q,pt,m,r,t,dra,dta,ur,T,ν,μ,eos; delta_f)*deg*fact
end

#spectra in a single pt point
function spectra_analytic(pt::C,fo::FreezeOutResult{A,B},part::particle_attribute{S,R,U,V}) where {C<:Number,A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)

    quadgk(alpha->_pointwise_spectra_analytic(pt,alpha,x,phi,part,eos),lb...,rb...)

end

#spectra points between max, min with given step (even spacing)
function spectra_analytic(fo::FreezeOutResult{A,B},part::particle_attribute{S,R,U,V},eos;pt_min=0.,pt_max=10.0,step=100,delta_f=false) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)

    #rb=min.(rb,rightbound)
    buff=alloc_segbuf(Float64, eltype(lb),Float64 ;size=1)
    
    [quadgk(alpha->_pointwise_spectra_analytic(pt,alpha,x,phi,part,eos;delta_f),lb...,rb...;segbuf=buff) for pt in range(pt_min,pt_max,step) ] 

end

#spectra points in a given range (can also be uneven spacing)
function spectra_analytic(fo::FreezeOutResult{A,B},part::particle_attribute{S,R,U,V}, pt_range::Vector{Float64};rtol=1000*sqrt(eps())) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)
    rb=min.(rb,rightbound)
    buff=alloc_segbuf(Float64, eltype(lb),Float64 ;size=1)
    [quadgk(alpha->_pointwise_spectra_analytic(pt,alpha,x,phi,part,eos),lb...,rb...;segbuf=buff,rtol=rtol) for pt in pt_range ] 

end


function spectra_analytic(fo::FreezeOutResult{A,B},part::particle_attribute{S,R,U,V},pt::C) where {C<:AbstractVector,A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)
    buff=alloc_segbuf(Float64, eltype(lb),Float64;size=1)

    [quadgk(alpha->_pointwise_spectra_analytic(i,alpha,x,phi,part,eos),lb...,rb...,segbuf=buff) for i in pt ] 

end


function multiplicity_analytic(pt,spectra)
    f = [spectra[i][1]*pt[i] for i in 1:length(pt)]
    2.0*π*NumericalIntegration.integrate(pt,f)
end

function multiplicity_analytic(fo::FreezeOutResult{A,B},part::particle_attribute{S,R,U,V},eos;rtol=1000*sqrt(eps()),delta_f = false,pt_min=0.,pt_max=10.0) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)
    #rb = min.(rb,(rightbound,))
    hcubature( b ->2.0*π *_pointwise_spectra_analytic(b[1],b[2],x,phi,part,eos;delta_f)*b[1],(pt_min,lb...),(pt_max,rb...);rtol=rtol)    
end

