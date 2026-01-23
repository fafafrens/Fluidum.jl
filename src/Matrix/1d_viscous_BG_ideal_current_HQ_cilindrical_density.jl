# the convention here are T, ur,  \[Pi]phiphi, \[Pi]etaeta, \[Pi]B, α, nur

function matrix1d_visc_BG_ideal_current_density_HQ!(A_i,Source,ϕ,t,X,params;free=true)
    
    #@show t, ϕ
    dpt = pressure_derivative(ϕ[1],Val(1),params.eos) #entropy
    dptt = pressure_derivative(ϕ[1],Val(2),params.eos)
        
    etaVisc = viscosity(ϕ[1], dpt, params.shear)
    tauShear = τ_shear(ϕ[1], dpt, params.shear)
    tauBulk = τ_bulk(ϕ[1], dpt, dptt, params.bulk)
    zetaBulkVisc = bulk_viscosity(ϕ[1], dpt, params.bulk)
    

    if free == true 
        # n is field number 6 here !
        Ds = diffusion(ϕ[1],ϕ[6],params.diffusion)
        tauDiff= τ_diffusion(ϕ[1],params.diffusion)
    else 
        #Ds = diffusion_hadron(ϕ[1],ϕ[6],params.eos,params.diffusion) #diffusion coefficient for hadrons
        #tauDiff=τ_diffusion_hadron(ϕ[1],ϕ[6],params.eos,params.diffusion) #tau diffusion for hadrons
    end

    
    dmp   = 0 #for now we don t have chemical potential in the eos
    dtdmp = 0 
    dmdmp = 0

    #actually our equations don t depend on p: we can just put as entry dpt instead, in any case it will not be used (but in the future maybe it will be )

    #(At,Ax, source)=one_d_viscous_matrix(ϕ,t,X[1],dpt,dpt,dptt,dmp,dtdmp,dmdmp,zeta,etaVisc,tauS,tauB,n,dtn,dmn,tauDiff,Ds)
    #(At,Ax, source)=one_d_viscous_matrix(ϕ,t,X[1],dpt,dpt,dptt,zeta,etaVisc,tauS,tauB,n,dtn,dmn,tauDiff,κ)
    #everything that is dpt where dpt not used it just to pass something
    (At,Ax, source)=one_d_viscous_matrix4(ϕ, t, X[1], dpt, dpt, dptt, zetaBulkVisc, etaVisc, tauShear, tauBulk, dpt, dpt, dpt, tauDiff, Ds)
        
    Ainv= inv(At)


    A_mul_B!(A_i[1], Ainv,Ax)
    
    
    jgemvavx!(Source, Ainv,source)
    
    end 

    
    function one_d_viscous_matrix4(X, tau, r, P_T, dP_dT, dP_dTdT, zetaBulkVisc, etaShearVisc, tauShear, tauBulk, n, dn_dT, dn_dmu, tauDiffusion, Ds)
    
        At=SMatrix{6,6}((X[2] .^2 .*dP_dT  .+ X[1] .*(1  .+ X[2] .^2) .*dP_dTdT,

X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ X[1] .*dP_dTdT),

0,

0,

0,

0,

2 .*X[2] .*( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT),

((1  .+ 2 .*X[2] .^2) .*( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT)) ./sqrt.(1  .+ X[2] .^2),

( .-2 .*etaShearVisc .*X[2]) ./(3. .*r .^2 .*sqrt.(1  .+ X[2] .^2)),

( .-2 .*etaShearVisc .*X[2]) ./(3. .*sqrt.(1  .+ X[2] .^2) .*tau .^2),

(zetaBulkVisc .*X[2]) ./sqrt.(1  .+ X[2] .^2),

0,

 .-X[2] .^2,

 .-(X[2] .*sqrt.(1  .+ X[2] .^2)),

(tauShear .*sqrt.(1  .+ X[2] .^2)) ./r .^2,

0,

0,

0,

 .-X[2] .^2,

 .-(X[2] .*sqrt.(1  .+ X[2] .^2)),

0,

(tauShear .*sqrt.(1  .+ X[2] .^2)) ./tau .^2,

0,

0,

X[2] .^2,

X[2] .*sqrt.(1  .+ X[2] .^2),

0,

0,

tauBulk .*sqrt.(1  .+ X[2] .^2),

0,

0,

0,

0,

0,

0,

sqrt.(1  .+ X[2] .^2))
        )
            #########################################################################################################################################################################
            
        Ax=SMatrix{6,6}((X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ X[1] .*dP_dTdT),

(1  .+ X[2] .^2) .*dP_dT  .+ X[1] .*X[2] .^2 .*dP_dTdT,

0,

0,

0,

( .-3 .*tauShear .*tauBulk .*X[1] .*X[2] .*X[6] .*dP_dT .*dP_dTdT) ./( .-3 .*tauShear .*tauBulk .*X[1] .*X[2] .^2 .*dP_dT .^2  .- X[1] .*(3 .*zetaBulkVisc .*tauShear .*X[2] .^2  .+ 4 .*etaShearVisc .*tauBulk .*X[2] .^2  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[3]  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[4]  .- 3 .*tauShear .*tauBulk .*X[5]  .- 3 .*tauShear .*tauBulk .*X[2] .^2 .*X[5]) .*dP_dTdT  .+ 3 .*tauShear .*tauBulk .*dP_dT .*(X[2] .^2 .*X[3]  .+ X[2] .^2 .*X[4]  .- X[2] .^2 .*X[5]  .+ X[1] .^2 .*dP_dTdT  .+ X[1] .^2 .*X[2] .^2 .*dP_dTdT)),

((1  .+ 2 .*X[2] .^2) .*( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT)) ./sqrt.(1  .+ X[2] .^2),

2 .*X[2] .*( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT),

( .-2 .*etaShearVisc) ./(3. .*r .^2),

( .-2 .*etaShearVisc) ./(3. .*tau .^2),

zetaBulkVisc,

(3 .*tauShear .*tauBulk .*X[1] .*X[6] .*( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT) .*dP_dTdT) ./( .-3 .*tauShear .*tauBulk .*X[1] .*X[2] .^2 .*dP_dT .^2  .- X[1] .*(3 .*zetaBulkVisc .*tauShear .*X[2] .^2  .+ 4 .*etaShearVisc .*tauBulk .*X[2] .^2  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[3]  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[4]  .- 3 .*tauShear .*tauBulk .*X[5]  .- 3 .*tauShear .*tauBulk .*X[2] .^2 .*X[5]) .*dP_dTdT  .+ 3 .*tauShear .*tauBulk .*dP_dT .*(X[2] .^2 .*X[3]  .+ X[2] .^2 .*X[4]  .- X[2] .^2 .*X[5]  .+ X[1] .^2 .*dP_dTdT  .+ X[1] .^2 .*X[2] .^2 .*dP_dTdT)),

 .-(X[2] .*sqrt.(1  .+ X[2] .^2)),

 .-1  .- X[2] .^2,

(tauShear .*X[2]) ./r .^2,

0,

0,

(3 .*tauShear .*tauBulk .*X[1] .*X[2] .*X[6] .*dP_dTdT) ./( .-3 .*tauShear .*tauBulk .*X[1] .*X[2] .^2 .*dP_dT .^2  .- X[1] .*(3 .*zetaBulkVisc .*tauShear .*X[2] .^2  .+ 4 .*etaShearVisc .*tauBulk .*X[2] .^2  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[3]  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[4]  .- 3 .*tauShear .*tauBulk .*X[5]  .- 3 .*tauShear .*tauBulk .*X[2] .^2 .*X[5]) .*dP_dTdT  .+ 3 .*tauShear .*tauBulk .*dP_dT .*(X[2] .^2 .*X[3]  .+ X[2] .^2 .*X[4]  .- X[2] .^2 .*X[5]  .+ X[1] .^2 .*dP_dTdT  .+ X[1] .^2 .*X[2] .^2 .*dP_dTdT)),

 .-(X[2] .*sqrt.(1  .+ X[2] .^2)),

 .-1  .- X[2] .^2,

0,

(tauShear .*X[2]) ./tau .^2,

0,

(3 .*tauShear .*tauBulk .*X[1] .*X[2] .*X[6] .*dP_dTdT) ./( .-3 .*tauShear .*tauBulk .*X[1] .*X[2] .^2 .*dP_dT .^2  .- X[1] .*(3 .*zetaBulkVisc .*tauShear .*X[2] .^2  .+ 4 .*etaShearVisc .*tauBulk .*X[2] .^2  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[3]  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[4]  .- 3 .*tauShear .*tauBulk .*X[5]  .- 3 .*tauShear .*tauBulk .*X[2] .^2 .*X[5]) .*dP_dTdT  .+ 3 .*tauShear .*tauBulk .*dP_dT .*(X[2] .^2 .*X[3]  .+ X[2] .^2 .*X[4]  .- X[2] .^2 .*X[5]  .+ X[1] .^2 .*dP_dTdT  .+ X[1] .^2 .*X[2] .^2 .*dP_dTdT)),

X[2] .*sqrt.(1  .+ X[2] .^2),

1  .+ X[2] .^2,

0,

0,

tauBulk .*X[2],

( .-3 .*tauShear .*tauBulk .*X[1] .*X[2] .*X[6] .*dP_dTdT) ./( .-3 .*tauShear .*tauBulk .*X[1] .*X[2] .^2 .*dP_dT .^2  .- X[1] .*(3 .*zetaBulkVisc .*tauShear .*X[2] .^2  .+ 4 .*etaShearVisc .*tauBulk .*X[2] .^2  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[3]  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[4]  .- 3 .*tauShear .*tauBulk .*X[5]  .- 3 .*tauShear .*tauBulk .*X[2] .^2 .*X[5]) .*dP_dTdT  .+ 3 .*tauShear .*tauBulk .*dP_dT .*(X[2] .^2 .*X[3]  .+ X[2] .^2 .*X[4]  .- X[2] .^2 .*X[5]  .+ X[1] .^2 .*dP_dTdT  .+ X[1] .^2 .*X[2] .^2 .*dP_dTdT)),

0,

0,

0,

0,

0,

X[2])
            )
            #########################################################################################################################################################################
        
        source=SVector{6}((( .-(X[2] .*X[3] .*(r .*X[2]  .+ sqrt.(1  .+ X[2] .^2) .*tau))  .+ X[4] .*(r  .- r .*X[2] .^2  .- X[2] .*sqrt.(1  .+ X[2] .^2) .*tau)  .+ (r  .+ r .*X[2] .^2  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau) .*(X[5]  .+ X[1] .*dP_dT)) ./(r .*tau),

( .-(X[4] .*(r .*X[2] .*sqrt.(1  .+ X[2] .^2)  .+ tau  .+ X[2] .^2 .*tau))  .- X[3] .*(r .*X[2] .*sqrt.(1  .+ X[2] .^2)  .+ 2 .*tau  .+ X[2] .^2 .*tau)  .+ X[2] .*(r .*sqrt.(1  .+ X[2] .^2)  .+ X[2] .*tau) .*(X[5]  .+ X[1] .*dP_dT)) ./(r .*tau),

( .-2 .*etaShearVisc .*r .*sqrt.(1  .+ X[2] .^2)  .+ 4 .*etaShearVisc .*X[2] .*tau  .+ 3 .*r .*X[3] .*tau) ./(3. .*r .^3 .*tau),

(4 .*etaShearVisc .*r .*sqrt.(1  .+ X[2] .^2)  .- 2 .*etaShearVisc .*X[2] .*tau  .+ 3 .*r .*X[4] .*tau) ./(3. .*r .*tau .^3),

(zetaBulkVisc .*X[2]) ./r  .+ X[5]  .+ (zetaBulkVisc .*sqrt.(1  .+ X[2] .^2)) ./tau,

( .-3 .*X[6] .*(X[1] .*(tauShear .*tauBulk .*r .*X[5]  .+ tauBulk .*r .*X[2] .^4 .*( .-2 .*etaShearVisc  .+ tauShear .*X[5])  .+ tauShear .*tauBulk .*X[2] .*sqrt.(1  .+ X[2] .^2) .*X[5] .*tau  .+ tauBulk .*X[2] .^3 .*sqrt.(1  .+ X[2] .^2) .*( .-2 .*etaShearVisc  .+ tauShear .*X[5]) .*tau  .+ r .*X[2] .^2 .*( .-2 .*etaShearVisc .*tauBulk  .+ tauShear .*X[5] .*(2 .*tauBulk  .+ sqrt.(1  .+ X[2] .^2) .*tau))  .+ tauShear .*tauBulk .*X[1] .*(1  .+ X[2] .^2) .*(r  .+ r .*X[2] .^2  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau) .*dP_dT) .*dP_dTdT  .+ tauBulk .*X[4] .*(tauShear .*X[2] .^2 .*(2 .*r .*(1  .+ X[2] .^2)  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau) .*dP_dT  .+ r .*X[1] .*(tauShear .*( .-1  .+ X[2] .^4)  .- X[2] .^2 .*sqrt.(1  .+ X[2] .^2) .*tau) .*dP_dTdT)  .+ tauBulk .*X[3] .*(tauShear .*X[2] .^2 .*(r  .+ r .*X[2] .^2  .+ 2 .*X[2] .*sqrt.(1  .+ X[2] .^2) .*tau) .*dP_dT  .- X[1] .*(r .*X[2] .^2 .*sqrt.(1  .+ X[2] .^2) .*tau  .+ tauShear .*(1  .+ X[2] .^2) .*(r  .- X[2] .*sqrt.(1  .+ X[2] .^2) .*tau)) .*dP_dTdT))) ./(r .*sqrt.(1  .+ X[2] .^2) .*tau .*(3 .*tauShear .*tauBulk .*X[2] .^2 .*dP_dT .*( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT)  .+ X[1] .*( .-3 .*tauShear .*tauBulk .*X[5]  .+ X[2] .^2 .*(3 .*zetaBulkVisc .*tauShear  .+ 4 .*etaShearVisc .*tauBulk  .- 3 .*tauShear .*tauBulk .*X[5])  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*(X[3]  .+ X[4]  .- X[1] .*dP_dT)) .*dP_dTdT)))
           )
            
        return (At,Ax, source)
        end
        
    