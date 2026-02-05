abstract type ParticleSpecies end

struct particle_simple{S} <:ParticleSpecies
    name::S
    mass::Float64
    degeneracy::Int64
    charge::Int64
    pt_list::Vector{Float64}
end

struct particle_full{T,S} <:ParticleSpecies
    name::S
    mass::Float64
    degeneracy::Int
    charge::Int
    pt_list::Vector{Float64}
    fj::T
end

name(particle::ParticleSpecies) = particle.name 
pt_list(particle::ParticleSpecies)= particle.pt_list 

pt_lists(species_list::Vector{T}) where {T<:ParticleSpecies} = vcat([pt_list(species) for species in species_list]...)


#particle_simple() = particle_simple("pion",0.13957,1,0,collect(0.1:0.5:2.0))
particle_simple(name, mass, degeneracy, charge) = particle_simple(name, mass, degeneracy, charge, collect(0.5:0.5:4.0))
particle_full(name, mass, degeneracy, charge,fj) = particle_full(name, mass, degeneracy, charge, collect(0.5:0.5:4.0),fj)


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
    return SVector(mt*cosh(eta_p-eta), pT*c, pT*s, zero(mt)) #last component is zero because it is always contracted with zero
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



@fastmath function dn_dpdx(fo,particle_species::particle_simple,coords,eta,pT, phi_p, eta_p, fun = exp)
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
    f_eq = fun(dot(pμ_up,uμ_down)/T)*exp(charge*α)
   
    deg*dot(dsigma_down(fo,coords),pμ_up)*f_eq/(2*pi)^3*Fluidum.fmGeV3
end



@fastmath function dn_dpdx(fo,particle_species::particle_full,coords,eta,pT, phi_p, eta_p)
    field = fo.fields
    fields_on_coords = field(coords...)
    @inbounds T = fields_on_coords[1]
    @inbounds ux = fields_on_coords[2]
    @inbounds uy = fields_on_coords[3]

    uμ_down = SVector(-sqrt(ux^2+uy^2+one(ux)),ux,uy,zero(T))
    uμ_up = SVector(sqrt(ux^2+uy^2+one(ux)),ux,uy,zero(T))

    m = particle_species.mass
    charge = particle_species.charge
    deg = particle_species.degeneracy
    f1 = particle_species.fj[1] #.. access Fj of particle
    f2 = particle_species.fj[2]
    #.. access Fj of particle

    @inbounds α = fields_on_coords[8]
    pμ_up = pmu_up(m, pT, phi_p, eta_p, eta)
    #f_eq = fun(dot(pμ_up,uμ_down)/T)*exp(charge*α)
    dσ_down = dsigma_down(fo,coords)
    Ep = -dot(pμ_up,uμ_down)
    su_contraction = Ep*dot(dσ_down,uμ_up)
    term1 = f1(Ep)*(dot(dσ_down,pμ_up)-su_contraction)
    term2 = f2(Ep)*su_contraction
    return (term1+term2)*exp(charge*α)/(2*pi)^3*Fluidum.fmGeV3
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


function dvn_dp_list_delta(fo,species_list, eta_p, wavenum_list; eta_min=-5.0, eta_max=5.0)

    pTlists = pt_list.(species_list)
    delta_lists = increment.(pTlists)
    domain = ([Fluidum.leftbounds(fo.x)...,eta_min,0.,-1.],[Fluidum.rightbounds(fo.x)...,eta_max,2pi,+1.])

    pt_length_max = maximum(length.(pt_list.(species_list)))
    function f(y,u,p) 
        fo, species_list, eta_p, wavenum_m = p
        @inbounds for k in eachindex(species_list)
            species = species_list[k]
            pTlist = pTlists[k]
            delta_list = delta_lists[k]
            @inbounds for i in eachindex(pTlist)
                delta = delta_list[i]
                pT = u[5]*delta + pTlist[i] #mapping from [0,1] to [pTlist[i]-delta, pTlist[i]+delta]
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
    prototype = zeros(3,pt_length_max,length(wavenum_list),length(species_list))
    par = (fo, species_list, eta_p, wavenum_list)
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




"""
q_vector_event_pt_dependent(result::ObservableResult,species_list, wavenum_list)


given a result::ObservableResult, computes the real q-vector as a function of pT for each species in species_list and wavenumber in wavenum_list
"""
function q_vector_event_pt_dependent(result::ObservableResult,species_list, wavenum_list)
    pTlists = pt_list.(species_list)
    pt_length_max = maximum(length.(pt_list.(species_list)))
    glauber, vn = result.glauber_multiplicity, result.vn
    qvec_result = zeros(pt_length_max,length(species_list),length(wavenum_list))
    for k in eachindex(species_list)
        pTlist = pTlists[k]
        for i in eachindex(pTlist)
            for wavenum in eachindex(wavenum_list)
                qvec_result[i,k,wavenum]+= sqrt(vn[1,i,wavenum,k]^2+vn[2,i,wavenum,k]^2)/vn[3,i,wavenum,k]
            end
        end
    end
    return qvec_result
end


function spectra_event(result::ObservableResult,species_list)
    pTlists = pt_list.(species_list)
    pt_length_max = maximum(length.(pt_list.(species_list)))
    glauber, vn = result.glauber_multiplicity, result.vn
    spectra_result = zeros(pt_length_max,length(species_list))
    for k in eachindex(species_list)
        pTlist = pTlists[k]
        for i in eachindex(pTlist)
            spectra_result[i,k] = vn[3,i,1,k]
        end
    end
    return spectra_result
end

function spectra(event_list,species_list)
    pTlists = pt_list.(species_list)
    pt_length_max = maximum(length.(pt_list.(species_list)))
    spectra_result = zeros(pt_length_max,length(species_list))
    for result in event_list
        glauber, vn = result.glauber_multiplicity, result.vn 
        for k in eachindex(species_list)
            pTlist = pTlists[k]
            for i in eachindex(pTlist)
                spectra_result[i,k] += vn[3,i,1,k]/length(event_list)
            end
        end
    end
    return spectra_result
end

"""
multiplicity_event(result::ObservableResult,species_list)


returns the total multiplicity M of all charged particles and identified particles in the species list for a given event result
"""
function multiplicity_event(result::ObservableResult,species_list)
    glauber, vn = result.glauber_multiplicity, result.vn
    pTlists = pt_list.(species_list)
    M = 0.
    M_species = zeros(length(species_list))
    for k in eachindex(species_list)
        pTlist = pTlists[k]
        for i in eachindex(pTlist)
            M+=vn[3,i,1,k] 
            M_species[k]+=vn[3,i,1,k]
        end
    end
    return (total_multiplicity = M, identified_multiplicity = M_species)
end

"""
g_species_event_pt_dependent(result::ObservableResult,species_list)


returns the multiplicity of each species in each pt bin normalized by the total multiplicity for a given event result
"""
function g_species_event_pt_dependent(result::ObservableResult,species_list)
    pTlists = pt_list.(species_list)
    pt_length_max = maximum(length.(pt_list.(species_list)))
    glauber, vn = result.glauber_multiplicity, result.vn
    g_result = zeros(pt_length_max,length(species_list))
    M_total = multiplicity_event(result,species_list).total_multiplicity
    for k in eachindex(species_list)
        pTlist = pTlists[k]
        for i in eachindex(pTlist)
            g_result[i,k]+= vn[3,i,1,k]/M_total
        end
    end
    return g_result  
end


"""
q_vector_event_integrated(result::ObservableResult,species_list, wavenum_list)


given a result::ObservableResult, computes the integrated q-vector for each wavenumber in wavenum_list
"""
function q_vector_event_integrated(result::ObservableResult,species_list, wavenum_list)
    pTlists = pt_list.(species_list)        
    glauber, vn = result.glauber_multiplicity, result.vn
    q_result = zeros(length(wavenum_list))
    g = g_species_event_pt_dependent(result,species_list)
    q_vector = q_vector_event_pt_dependent(result,species_list,wavenum_list)
    for k in eachindex(species_list)
        pTlist = pTlists[k]
        for i in eachindex(pTlist)
            q_result += g[i,k] * q_vector[i,k,:]
        end
    end
    return q_result
end

"""
harmonic_coefficient(event_list,species_list, wavenum_list)


computes the harmonic coefficients v_n{2} for a list of events, for each species in species_list and wavenumber in wavenum_list
"""
function harmonic_coefficient(event_list,species_list, wavenum_list)
    pTlists = pt_list.(species_list)        
    pt_length_max = maximum(length.(pt_list.(species_list)))
    
    vm_result = zeros(pt_length_max,length(species_list),length(wavenum_list))
    vm_result_integrated = zeros(length(species_list),length(wavenum_list))
    vm_result_charged = zeros(pt_length_max,length(wavenum_list))
    vm_result_charged_integrated = zeros(length(wavenum_list))
        
    for result in event_list        
        q_vector_pt = q_vector_event_pt_dependent(result,species_list,wavenum_list)
        q_vector_total = q_vector_event_integrated(result,species_list,wavenum_list)

        for k in eachindex(species_list)
            pTlist = pTlists[k]
            for i in eachindex(pTlist)
                for wavenum in eachindex(wavenum_list)
                    vm_result[i,k,wavenum] += q_vector_pt[i,k,wavenum]*q_vector_total[wavenum]/length(event_list)
                end
            end
        end

    end

    vm_result_integrated = sum(vm_result, dims=1)
    vm_result_charged = sum(vm_result, dims=2)
    vm_result_charged_integrated = sum(vm_result_charged, dims=1)

    return (vm_result=vm_result, vm_result_integrated=vm_result_integrated, vm_result_charged=vm_result_charged, vm_result_charged_integrated=vm_result_charged_integrated)

   
end


"""
select_cc_events_glauber(data,cc_fraction)

divides data into cc_fraction centrality classes based on glauber multiplicity
"""
function select_cc_events_glauber(data,cc_fraction)
    glauber_vec = extract_glauber_multiplicity(data)
    sorted_indices = sortperm(glauber_vec,rev=true)
    Nev = length(glauber_vec)
    cc_ev_num = div(cc_fraction*Nev,100)
    selected_data = [data[sorted_indices[1+i*cc_ev_num:cc_ev_num+i*cc_ev_num]] for i in 0:cc_fraction-1]
    return selected_data
end
