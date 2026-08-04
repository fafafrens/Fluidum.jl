@inline function _dn_pTdpTdx_internal_spectra(pT, ϕ_p, η_p, α, ϕ_s, η_s, x::A, fields::B, part::Particle{R,N,S}, eos, transport; delta_f = false) where {A<:SplineInterp,B<:SplineInterp, R, N, S}

    τ,r= x(α)
    dta,dra=jacobian(x,α)
    T,ur,pi_phi,pi_eta,pi_b,fugacity,νr=fields(α) 
    
    m = part.Mass
    charm = part.Nc+part.Nac
    degeneracy = part.Degeneracy

    if charm != 0
        statistics = 0 #boltzmann
    else
        if iseven(2*part.Spin)
            statistics = -1 #fermions
        else
            statistics = 1 #bosons
        end
    end

    fact = besseli(1, eos.hadron_list.ccbar/2)/besseli(0, eos.hadron_list.ccbar/2)

    if charm == 2
        fact = 1
    end 

    mT = hypot(pT,m)
    uτ = hypot(ur, 1)
    chη = cosh(η_p-η_s)
    cϕ = cos(ϕ_p-ϕ_s)
    pdotsigma = dra*mT*chη-dta*pT*cϕ
    det_g = -τ*r
    Ebar = mT*chη*uτ - pT*cϕ*ur
    
    charm_prefactor = -inv(T)*charm/normalization(T,fugacity,eos)
    bulk_prefactor = zero(charm_prefactor)
    shear_prefactor = zero(charm_prefactor)
    entropy = pressure_derivative(T,Val(1),eos)
    dtdtp = pressure_derivative(T, Val(2), eos)
    bulk, shear = transport
    cs2 = entropy / (T * dtdtp)

    if iszero(charm)
        bulk_prefactor = τ_bulk(T,entropy,dtdtp,bulk)/bulk_viscosity(T,entropy,bulk)*(Ebar/T*(inv(3.)-cs2)-inv(3.)*m^2/T/Ebar)
        shear_prefactor = 0.5/T^3/entropy
    end

    f0 = inv(exp(+uτ*mT/T*chη-ur*pT/T*cϕ)+statistics)*exp(charm*fugacity)
    
    if delta_f
        δf = (-mT*chη*νr*ur/uτ+pT*cϕ*νr)*charm_prefactor +
             (one(f0)+statistics*f0)*bulk_prefactor*pi_b +
             (one(f0)+statistics*f0)*shear_prefactor*(-(pi_phi+pi_eta)*(mT*cosh(η_p)*ur-pT*cos(ϕ_p)*uτ)^2+pi_phi*(-pT*sin(ϕ_p))^2+pi_eta*(mT^2*sinh(η_p))^2)
    else
        δf = zero(f0)
    end

    f = f0*(one(f0)+ δf)

    return (2pi)/(2pi)^3*degeneracy*f*det_g*pdotsigma*fmGeV^3
end

#spectra in a single pt point
#numerical integration of function of alpha, from the lb of alpha to the rb of alpha (over all the fo surface)
function spectra(pT::C,fo::FreezeOutResult{A,B},part::Particle{R,N,S},eos, transport;rtol=1000*sqrt(eps()), η_s_max = 5., delta_f = false ) where {C<:Number,A<:SplineInterp,B<:SplineInterp,R,N,S}
    x,fields = fo   
    lb=leftbounds(x)
    rb=rightbounds(x)
    η_p = 0.
    domain = ([lb...,0.,0.,0,],[rb...,η_s_max,2π,2π])
    function f(u,p)
        fo, part, η_p = p
        2. *_dn_pTdpTdx_internal_spectra(pT, u[3], η_p, u[1], u[4], u[2], x, fields, part, eos, transport; delta_f = delta_f)
    end
    prototype = zero(pT)
    par = (fo, part, η_p)
    prob = IntegralProblem(IntegralFunction(f,prototype),domain,par)
    result = solve(prob, CubaVegas(), reltol=rtol, abstol=1e-6)
    return result
end
