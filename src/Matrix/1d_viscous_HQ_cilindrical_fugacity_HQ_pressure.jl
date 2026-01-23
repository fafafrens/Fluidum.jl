# the convention here are T, ur,  \[Pi]phiphi, \[Pi]etaeta, \[Pi]B, α, nur

function matrix1d_visc_HQ_pressure!(A_i,Source,ϕ,t,X,params;free=true)
    
    
    pressure_field = hq_pressure(ϕ[1],ϕ[6];m=params.diffusion.mass)
    dPhq_dT, dPhq_dmu = pressure_field.pressure_derivative
    dPhq_dTdT,dPhq_dTdmu,ddp_dmudmu_HQ = pressure_field.pressure_hessian

    #@show t, ϕ
    dpt = pressure_derivative(ϕ[1], Val(1), params.eos) #entropy
    dptt = pressure_derivative(ϕ[1], Val(2), params.eos)

    dpt_tot = dpt + dPhq_dT
    dptt_tot = dptt + dPhq_dTdT

    #dPhq_dT += 1e-5
    #dPhq_dmu += 1e-5
    #dPhq_dTdT += 1e-5
    #dPhq_dTdmu += 1e-5
    #ddp_dmudmu_HQ += 1e-5
    #dpt_tot += 1e-10
    #dptt_tot += 1e-10
    etaVisc = viscosity(ϕ[1], dpt_tot, params.shear)
    tauShear = τ_shear(ϕ[1], dpt_tot, params.shear)
    tauB    = τ_bulk(ϕ[1], dpt_tot, dptt_tot, params.bulk)
    zeta    = bulk_viscosity(ϕ[1], dpt_tot, params.bulk)
    if free == true 
        thermo = hq_density(ϕ[1],ϕ[6];m=params.diffusion.mass)
            n = thermo.value
            dtn, dmn = thermo.gradient
    else 
        thermo = thermodynamic(ϕ[1],ϕ[6],params.eos.hadron_list)
            n = thermo.pressure
    dtn, dmn = thermo.pressure_derivative
    end
    

    dmn_eps = let s = get(ENV, "FLUIDUM_DMN_EPS", "1e-4")
        v = tryparse(Float64, s)
        v === nothing ? 1e-4 : v
    end
    dtn_eps = let s = get(ENV, "FLUIDUM_DTN_EPS", "1e-4")
        v = tryparse(Float64, s)
        v === nothing ? 1e-4 : v
    end

    dmn += dmn_eps
    dtn += dtn_eps

    
    if free == true 
        Ds = diffusion(ϕ[1],n,params.diffusion)
        tauDiff= τ_diffusion(ϕ[1],params.diffusion)
    else 
        Ds = diffusion_hadron(ϕ[1],ϕ[6],params.eos,params.diffusion) #diffusion coefficient for hadrons
        tauDiff=τ_diffusion_hadron(ϕ[1],ϕ[6],params.eos,params.diffusion) #tau diffusion for hadrons
    end

    (At,Ax, source)=one_d_viscous_matrix13(ϕ,t,X[1],nothing,dpt,dptt,dPhq_dT,dPhq_dmu, dPhq_dTdT,dPhq_dTdmu,ddp_dmudmu_HQ,zeta,etaVisc,tauShear,tauB,n,dtn,dmn,tauDiff,Ds)
        
    Ainv= inv(At)
    
    A_mul_B!(A_i[1], Ainv,Ax)
    
    jgemvavx!(Source, Ainv,source)
    
end 

    
    function one_d_viscous_matrix13(X,tau,r,pt,dP_dT,dP_dTdT,dPhq_dT,dPhq_dmu, dPhq_dTdT,dPhq_dTdmu,dPhq_dmudmu,zetaBulkVisc,etaShearVisc, tauShear, tauBulk, n, dn_dT, dn_dmu, tauDiffusion, Ds)
    
        At=SMatrix{7,7}(
     (X[2] .^2 .*dP_dT  .+ X[2] .^2 .*dPhq_dT  .+ (1  .+ X[2] .^2) .*(X[1] .*dP_dTdT  .+ X[6] .*dPhq_dTdmu  .+ X[1] .*dPhq_dTdT),

X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ X[1] .*dP_dTdT  .+ dPhq_dT  .+ X[6] .*dPhq_dTdmu  .+ X[1] .*dPhq_dTdT),

0,

0,

0,

sqrt.(1  .+ X[2] .^2) .*dn_dT,

0,

2 .*X[2] .*( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT  .+ X[6] .*dPhq_dmu  .+ X[1] .*dPhq_dT),

((1  .+ 2 .*X[2] .^2) .*( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT  .+ X[6] .*dPhq_dmu  .+ X[1] .*dPhq_dT)) ./sqrt.(1  .+ X[2] .^2),

( .-2 .*etaShearVisc .*X[2]) ./(3. .*r .^2 .*sqrt.(1  .+ X[2] .^2)),

( .-2 .*etaShearVisc .*X[2]) ./(3. .*sqrt.(1  .+ X[2] .^2) .*tau .^2),

(zetaBulkVisc .*X[2]) ./sqrt.(1  .+ X[2] .^2),

(n .*(X[2]  .+ X[2] .^3)  .+ X[7]) ./(1  .+ X[2] .^2) .^1.5,

 .-((tauDiffusion .*X[2] .*X[7]) ./sqrt.(1  .+ X[2] .^2)),

 .-X[2] .^2,

 .-(X[2] .*sqrt.(1  .+ X[2] .^2)),

(tauShear .*sqrt.(1  .+ X[2] .^2)) ./r .^2,

0,

0,

0,

0,

 .-X[2] .^2,

 .-(X[2] .*sqrt.(1  .+ X[2] .^2)),

0,

(tauShear .*sqrt.(1  .+ X[2] .^2)) ./tau .^2,

0,

0,

0,

X[2] .^2,

X[2] .*sqrt.(1  .+ X[2] .^2),

0,

0,

tauBulk .*sqrt.(1  .+ X[2] .^2),

0,

0,

X[2] .^2 .*dPhq_dmu  .+ (1  .+ X[2] .^2) .*(X[6] .*dPhq_dmudmu  .+ X[1] .*dPhq_dTdmu),

X[2] .*sqrt.(1  .+ X[2] .^2) .*(dPhq_dmu  .+ X[6] .*dPhq_dmudmu  .+ X[1] .*dPhq_dTdmu),

0,

0,

0,

sqrt.(1  .+ X[2] .^2) .*dn_dmu,

Ds .*n .*X[2] .*sqrt.(1  .+ X[2] .^2),

0,

0,

0,

0,

0,

X[2] ./sqrt.(1  .+ X[2] .^2),

tauDiffusion .*sqrt.(1  .+ X[2] .^2))
        )
            #########################################################################################################################################################################
            
        Ax=SMatrix{7,7}(
     (X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ X[1] .*dP_dTdT  .+ dPhq_dT  .+ X[6] .*dPhq_dTdmu  .+ X[1] .*dPhq_dTdT),

(1  .+ X[2] .^2) .*dP_dT  .+ (1  .+ X[2] .^2) .*dPhq_dT  .+ X[2] .^2 .*(X[1] .*dP_dTdT  .+ X[6] .*dPhq_dTdmu  .+ X[1] .*dPhq_dTdT),

0,

0,

0,

X[2] .*dn_dT,

0,

((1  .+ 2 .*X[2] .^2) .*( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT  .+ X[6] .*dPhq_dmu  .+ X[1] .*dPhq_dT)) ./sqrt.(1  .+ X[2] .^2),

2 .*X[2] .*( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT  .+ X[6] .*dPhq_dmu  .+ X[1] .*dPhq_dT),

( .-2 .*etaShearVisc) ./(3. .*r .^2),

( .-2 .*etaShearVisc) ./(3. .*tau .^2),

zetaBulkVisc,

n,

 .-((tauDiffusion .*X[2] .^2 .*X[7]) ./(1  .+ X[2] .^2)),

 .-(X[2] .*sqrt.(1  .+ X[2] .^2)),

 .-1  .- X[2] .^2,

(tauShear .*X[2]) ./r .^2,

0,

0,

0,

0,

 .-(X[2] .*sqrt.(1  .+ X[2] .^2)),

 .-1  .- X[2] .^2,

0,

(tauShear .*X[2]) ./tau .^2,

0,

0,

0,

X[2] .*sqrt.(1  .+ X[2] .^2),

1  .+ X[2] .^2,

0,

0,

tauBulk .*X[2],

0,

0,

X[2] .*sqrt.(1  .+ X[2] .^2) .*(dPhq_dmu  .+ X[6] .*dPhq_dmudmu  .+ X[1] .*dPhq_dTdmu),

(1  .+ X[2] .^2) .*dPhq_dmu  .+ X[2] .^2 .*(X[6] .*dPhq_dmudmu  .+ X[1] .*dPhq_dTdmu),

0,

0,

0,

X[2] .*dn_dmu,

Ds .*n .*(1  .+ X[2] .^2),

0,

0,

0,

0,

0,

1,

tauDiffusion .*X[2])
            )
            #########################################################################################################################################################################
        
        source=SVector{7}(
(( .-(X[2] .*X[3] .*(r .*X[2]  .+ sqrt.(1  .+ X[2] .^2) .*tau))  .+ X[4] .*(r  .- r .*X[2] .^2  .- X[2] .*sqrt.(1  .+ X[2] .^2) .*tau)  .+ (r  .+ r .*X[2] .^2  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau) .*(X[5]  .+ X[6] .*dPhq_dmu  .+ X[1] .*(dP_dT  .+ dPhq_dT))) ./(r .*tau),

( .-(X[4] .*(r .*X[2] .*sqrt.(1  .+ X[2] .^2)  .+ tau  .+ X[2] .^2 .*tau))  .- X[3] .*(r .*X[2] .*sqrt.(1  .+ X[2] .^2)  .+ 2 .*tau  .+ X[2] .^2 .*tau)  .+ X[2] .*(r .*sqrt.(1  .+ X[2] .^2)  .+ X[2] .*tau) .*(X[5]  .+ X[6] .*dPhq_dmu  .+ X[1] .*(dP_dT  .+ dPhq_dT))) ./(r .*tau),

( .-2 .*etaShearVisc .*r .*sqrt.(1  .+ X[2] .^2)  .+ 4 .*etaShearVisc .*X[2] .*tau  .+ 3 .*r .*X[3] .*tau) ./(3. .*r .^3 .*tau),

(4 .*etaShearVisc .*r .*sqrt.(1  .+ X[2] .^2)  .- 2 .*etaShearVisc .*X[2] .*tau  .+ 3 .*r .*X[4] .*tau) ./(3. .*r .*tau .^3),

(zetaBulkVisc .*X[2]) ./r  .+ X[5]  .+ (zetaBulkVisc .*sqrt.(1  .+ X[2] .^2)) ./tau,

X[7] .*(1 ./r  .+ X[2] ./(sqrt.(1  .+ X[2] .^2) .*tau))  .+ n .*(X[2] ./r  .+ sqrt.(1  .+ X[2] .^2) ./tau),

X[7])

        )
            
        return (At,Ax, source)
        end
        
