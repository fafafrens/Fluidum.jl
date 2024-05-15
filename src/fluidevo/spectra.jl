
function cilindrical_thermal_spectra(pt,m,r,t,dra,dta,ur,T)
    fmGeV = 5.068
    mt=sqrt(m^2+pt^2)
    factor=1/(2*pi^2)*t*r
    karg=mt*sqrt(1+ur^2)/T
    iarg=pt*ur/T

    r_factor=-factor* dra*mt
    t_factor=factor* dta*pt

    i0=besseli0(iarg)
    i1=besseli1(iarg)
    k0=besselk0(karg)
    k1=besselk1(karg)

    result= r_factor*k1*i0+t_factor*k0*i1
    return result*fmGeV^3
end 



function _pointwise_spectra(pt,m,alpha,x::A,phi::B) where {A<:SplineInterp,B<:SplineInterp}
    t,r= x(alpha)
    dta,dra=jacobian(x,alpha)
    T,ur,pi_phi,pi_eta,pi_b=phi(alpha)
    cilindrical_thermal_spectra(pt,m,r,t,dra,dta,ur,T)
end

#spectra in a single pt point
function spectra_analitic(pt::C,m::D,fo::FreezeOutResult{A,B}) where {C<:Number,D<:Number,A<:SplineInterp,B<:SplineInterp}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)

    #numerical integration of function of alpha, from the lb of alpha to the rb of alpha (over all the fo surface)
    #returns a pair (integral, error)
    quadgk(alpha->_pointwise_spectra(pt,m,alpha,x,phi),lb...,rb...)

end

#spectra points over pt range
function spectra_analitic(m::Number,fo::FreezeOutResult{A,B};pt_min=0.,pt_max=8.0,step=100) where {A<:SplineInterp,B<:SplineInterp}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)
    buff=alloc_segbuf(Float64, eltype(lb),Float64 ;size=1)
    
    [quadgk(alpha->_pointwise_spectra(pt,m,alpha,x,phi),lb...,rb...;segbuf=buff) for pt in range(pt_min,pt_max,step) ] 

end

#spectra points in a given range (can also be uneven spacing)
function spectra_analitic(fo::FreezeOutResult{A,B},part::particle_attribute{S,R,U,V}, pt_range::Vector{Float64};rtol=1000*sqrt(eps()),decays=true,rightbound=(100,),ccbar = 2.76) where {A<:SplineInterp,B<:SplineInterp,S,R,U,V}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)
    rb=min.(rb,rightbound)
    buff=alloc_segbuf(Float64, eltype(lb),Float64 ;size=1)
    m = part.mass
    [quadgk(alpha->_pointwise_spectra(pt,m,alpha,x,phi),lb...,rb...;segbuf=buff,rtol=rtol) for pt in pt_range ] 

end


#forgot mass?
function spectra_analitic(pt::C,m::Number,fo::FreezeOutResult{A,B}) where {C<:AbstractVector,A<:SplineInterp,B<:SplineInterp}
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)
    buff=alloc_segbuf(Float64, eltype(lb),Float64;size=1)

    [quadgk(alpha->_pointwise_spectra(i,m,alpha,x,phi),lb...,rb...,segbuf=buff) for i in pt ] 

end


function multiplicity_analitic(pt,spectra)
    f = [spectra[i][1]*pt[i] for i in 1:length(pt)]
    2.0*π*NumericalIntegration.integrate(pt,f)
end

function multiplicity_analitic(fo::FreezeOutResult{A,B},part::particle_attribute{S,R,U,V},kernel::kernels{T};rtol=1000*sqrt(eps())) where {A<:SplineInterp,B<:SplineInterp,T,S,R,U,V}
    m = part.mass
    x,phi=fo
    lb=leftbounds(x)
    rb=rightbounds(x)

    hcubature( b ->2.0*π *_pointwise_spectra(b[1],m,b[2],x,phi,kernel)*b[1],(0,lb...),(3,rb...);rtol=rtol)
    
end

