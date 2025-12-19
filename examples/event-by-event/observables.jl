using StaticArrays
using LinearAlgebra
using Integrals
using Cuba

struct particle_simple
    mass::Float64
    degeneracy::Int64
    charge::Int64
end

@inline @inbounds @fastmath function dsigma_down(fo, coords)
    point= fo.x
    t,x,y= point(coords...)
    jmatrix = Fluidum.jacobian(point,coords)
    
    
    @muladd a0=-t*jmatrix[1,2]*jmatrix[2,3]+t*jmatrix[1,3]*jmatrix[2,2]
    @muladd a1=-t*jmatrix[1,3]*jmatrix[2,1]+t*jmatrix[1,1]*jmatrix[2,3]
    @muladd a2=-t*jmatrix[1,1]*jmatrix[2,2]+t*jmatrix[1,2]*jmatrix[2,1]
    a3=zero(a0)

    return SVector(a0,a1,a2,a3) #with minus!!! #with determinant of the metric
    
end

@inline @fastmath function pmu_up(m, pT, phi_p, eta_p, eta)
    mt = sqrt(m^2 + pT^2)
    s,c=sincos(phi_p)
    return SVector(mt*cosh(eta_p-eta), pT*c, pT*s, zero(mt))
end



@inbounds function dn_dpdx(fo,m,coords,eta,pT, phi_p, eta_p)
    field = fo.fields
    fields_on_coords = field(coords...)
    T = fields_on_coords[1]
    ux = fields_on_coords[2]
    uy = fields_on_coords[3]

    uμ_down = SVector{4}(-sqrt(ux^2+uy^2+1),ux,uy,zero(T))
    α = fields_on_coords[8]
    pμ_up = pmu_up(m, pT, phi_p, eta_p, eta)
    f_eq = exp(dot(pμ_up,uμ_down)/T + α)
   
    dot(dsigma_down(fo,coords),pmu_up(m,pT,phi_p,eta_p,eta))*f_eq/(2*pi)^3*Fluidum.fmGeV^3
end



@fastmath function dn_dpdx(fo,particle_species::particle_simple,coords,eta,pT, phi_p, eta_p)
    field = fo.fields
    fields_on_coords = field(coords...)
    @inbounds T = fields_on_coords[1]
    @inbounds ux = fields_on_coords[2]
    @inbounds uy = fields_on_coords[3]

    uμ_down = SVector(-sqrt(ux^2+uy^2+one(ux)),ux,uy,zero(T))

    m = particle_species.mass
    charge = particle_species.charge
    deg = particle_species.degeneracy

    @inbounds α = fields_on_coords[8]
    pμ_up = pmu_up(m, pT, phi_p, eta_p, eta)
    f_eq = exp(dot(pμ_up,uμ_down)/T + charge*α)
   
    deg*dot(dsigma_down(fo,coords),pμ_up)*f_eq/(2*pi)^3*Fluidum.fmGeV3
end


function dn_dp(fo,m,pT, phi_p, eta_p; eta_min=-5.0, eta_max=5.0)
    domain = ([Fluidum.leftbounds(fo.x)...,eta_min],[Fluidum.rightbounds(fo.x)...,eta_max])
    function f(u,p) 
        fo, m, pT, phi_p, eta_p = p
        dn_dpdx(fo,m,(u[1],u[2]),u[3],pT, phi_p, eta_p)
    end
    p = (fo, m, pT, phi_p, eta_p)
    prob = IntegralProblem(f,domain,p)
    result = solve(prob, HCubatureJL(), reltol=1e-3, abstol=1e-6)
    return result
end


function dn_dpTdetap(fo,m,pT::Float64, eta_p; eta_min=-5.0, eta_max=5.0)
    domain = ([Fluidum.leftbounds(fo.x)...,eta_min,0],[Fluidum.rightbounds(fo.x)...,eta_max,2pi])
    function f(u,p) 
        fo, m, pT, eta_p = p
        dn_dpdx(fo,m,(u[1],u[2]),u[3],pT, u[4], eta_p)
    end
    p = (fo, m, pT, eta_p)
    prob = IntegralProblem(f,domain,p)
    result = solve(prob, HCubatureJL(), reltol=1e-3, abstol=1e-6)
    return result
end

function dn_dpTdetap(fo,m,pTlist, eta_p; eta_min=-5.0, eta_max=5.0)
    domain = ([Fluidum.leftbounds(fo.x)...,eta_min,0],[Fluidum.rightbounds(fo.x)...,eta_max,2pi])
    function f(y,u,p) 
        fo, m, pTlist, eta_p = p
        for (i,pT) in enumerate(pTlist)
        y[i]=dn_dpdx(fo,m,(u[1],u[2]),u[3],pT, u[4], eta_p)
        end
    end
    prototype = zeros(length(pTlist))
    p = (fo, m, pTlist, eta_p)
    prob = IntegralProblem(IntegralFunction(f,prototype),domain,p)
    result = solve(prob, HCubatureJL(), reltol=1e-3, abstol=1e-6)
    return result
end

function dvn_dp_list(fo,m, pTlist, eta_p, wavenum_list; eta_min=-5.0, eta_max=5.0)
    domain = ([Fluidum.leftbounds(fo.x)...,eta_min,0.],[Fluidum.rightbounds(fo.x)...,eta_max,2pi])
    function f(y,u,p) 
        fo, m, pTlist, eta_p, wavenum_m = p
    
        for (i,pT) in enumerate(pTlist)
            denom=dn_dpdx(fo,m,(u[1],u[2]),u[3],pT, u[4], eta_p)*pT #should this be the integrated yield or not? if yes, pt = u[5]
       
            for (j,wavenum_m) in enumerate(wavenum_list)
        y[1,i,j]=dn_dpdx(fo,m,(u[1],u[2]),u[3],pT, u[4], eta_p)*cos(wavenum_m*u[4])*pT
        y[2,i,j]=dn_dpdx(fo,m,(u[1],u[2]),u[3],pT, u[4], eta_p)*sin(wavenum_m*u[4])*pT
        y[3,i,j]=denom
        end
    end
    end
    prototype = zeros(3,length(pTlist),length(wavenum_list))
    par = (fo, m, pTlist, eta_p, wavenum_list)
    prob = IntegralProblem(IntegralFunction(f,prototype),domain,par)
    result = solve(prob, CubaVegas(), reltol=1e-3, abstol=1e-6)
    return result
end

function increment(list)
       de=diff(list)/2
       push!(de,last(de))
       return de
end

function indicator(x,qT,delta)
    if abs(x - qT) <= delta/2
        return 1.0
    else
        return 0.0
    end
end



function dvn_dp_list_delta(fo,species_list, pTlist, eta_p, wavenum_list; eta_min=-5.0, eta_max=5.0)

    delta_list = increment(pTlist)
    domain = ([Fluidum.leftbounds(fo.x)...,eta_min,0.,-1.],[Fluidum.rightbounds(fo.x)...,eta_max,2pi,+1.])

    function f(y,u,p) 
        fo, species_list, pTlist, eta_p, wavenum_m = p
        @inbounds for k in eachindex(species_list)
            species = species_list[k]
            @inbounds for i in eachindex(pTlist)
                delta = delta_list[i]
                pT = u[5]*delta + pTlist[i] #mapping from [0,1] to [pTlist[i], pTlist[i]+delta]
                denom=dn_dpdx(fo,species,(u[1],u[2]),u[3],pT, u[4], eta_p)*pT 
                @inbounds for j in eachindex(wavenum_list)
                    wavenum_m=wavenum_list[j]
                    snphi,cnphi=sincos(wavenum_m*u[4])
                    y[1,i,j,k]=denom*cnphi
                    y[2,i,j,k]=denom*snphi
                    y[3,i,j,k]=denom
                end
            end
        end
    end
    prototype = zeros(3,length(pTlist),length(wavenum_list),length(species_list))
    par = (fo, species_list, pTlist, eta_p, wavenum_list)
    prob = IntegralProblem(IntegralFunction(f,prototype),domain,par)
    result = solve(prob, CubaVegas(), reltol=1e-3, abstol=1e-6)
    return result
end


function dn_detap(fo,m,eta_p; eta_min=-5.0, eta_max=5.0,pt_min=0., pt_max=5.0)
    domain = ([Fluidum.leftbounds(fo.x)...,eta_min,0,pt_min],[Fluidum.rightbounds(fo.x)...,eta_max,2pi,pt_max])
    function f(u,p) 
        fo, m, eta_p = p
        dn_dpdx(fo,m,(u[1],u[2]),u[3],u[5], u[4], eta_p) * u[5]
    end
    p = (fo, m, eta_p)
    prob = IntegralProblem(f,domain,p)
    result = solve(prob, HCubatureJL(), reltol=1e-3, abstol=1e-6)
    return result
end
#=
function dvn_dp_flow(fo,m,pT, eta_p, wavenum_m; eta_min=-5.0, eta_max=5.0)
    domain = ([Fluidum.leftbounds(fo.x)...,eta_min,0.],[Fluidum.rightbounds(fo.x)...,eta_max,2pi])
    function f(y,u,p) 
        fo, m, pT, eta_p, wavenum_m = p
        y[1]=dn_dpdx(fo,m,(u[1],u[2]),u[3],pT, u[4], eta_p)*cos(wavenum_m*u[4])
        y[2]=dn_dpdx(fo,m,(u[1],u[2]),u[3],pT, u[4], eta_p)*sin(wavenum_m*u[4])
        y[3]=dn_dpdx(fo,m,(u[1],u[2]),u[3],pT, u[4], eta_p)
    end
    prototype = SVector{3}(0.,0.,0.)
    par = (fo, m, pT, eta_p, wavenum_m)
    prob = IntegralProblem(IntegralFunction(f,prototype),domain,par)
   # result = solve(prob, HCubatureJL(), reltol=1e-3, abstol=1e-6)
   # return result
end

=#

#bello! ptloop is inside
function dvn_dp_list(fo,m, pTlist, eta_p, wavenum_m; eta_min=-5.0, eta_max=5.0)
    domain = ([Fluidum.leftbounds(fo.x)...,eta_min,0.],[Fluidum.rightbounds(fo.x)...,eta_max,2pi])
    function f(y,u,p) 
        fo, m, pTlist, eta_p, wavenum_m = p
        for (i,pT) in enumerate(pTlist)
        y[1,i]=dn_dpdx(fo,m,(u[1],u[2]),u[3],pT, u[4], eta_p)*cos(wavenum_m*u[4])
        y[2,i]=dn_dpdx(fo,m,(u[1],u[2]),u[3],pT, u[4], eta_p)*sin(wavenum_m*u[4])
        y[3,i]=dn_dpdx(fo,m,(u[1],u[2]),u[3],pT, u[4], eta_p) #same as above. If integrated yield is needed, pt = u[5]
        end
    end
    prototype = zeros(3,length(pTlist))
    par = (fo, m, pTlist, eta_p, wavenum_m)
    prob = IntegralProblem(IntegralFunction(f,prototype),domain,par)
    result = solve(prob, CubaVegas(), reltol=1e-3, abstol=1e-6)
    return result
end




function dvn_dp_list(fo,m, pTlist, eta_p, wavenum_list; eta_min=-5.0, eta_max=5.0)
    domain = ([Fluidum.leftbounds(fo.x)...,eta_min,0.],[Fluidum.rightbounds(fo.x)...,eta_max,2pi])
    function f(y,u,p) 
        fo, m, pTlist, eta_p, wavenum_m = p
    
        for (i,pT) in enumerate(pTlist)
            denom=dn_dpdx(fo,m,(u[1],u[2]),u[3],pT, u[4], eta_p) #same as above. If integrated yield is needed, pt = u[5]
       
            for (j,wavenum_m) in enumerate(wavenum_list)
        y[1,i,j]=dn_dpdx(fo,m,(u[1],u[2]),u[3],pT, u[4], eta_p)*cos(wavenum_m*u[4])
        y[2,i,j]=dn_dpdx(fo,m,(u[1],u[2]),u[3],pT, u[4], eta_p)*sin(wavenum_m*u[4])
        y[3,i,j]=denom
        end
    end
    end
    prototype = zeros(3,length(pTlist),length(wavenum_list))
    par = (fo, m, pTlist, eta_p, wavenum_list)
    prob = IntegralProblem(IntegralFunction(f,prototype),domain,par)
    result = solve(prob, CubaVegas(), reltol=1e-3, abstol=1e-6)
    return result
end
