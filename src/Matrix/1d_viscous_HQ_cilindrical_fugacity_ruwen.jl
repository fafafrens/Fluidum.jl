# the convention here are T, ur,  \[Pi]phiphi, \[Pi]etaeta, \[Pi]B, α, nur

function matrix1d_visc_HQ_ruwen2!(A_i, Source, ϕ, t, X, params;free=true)

    #@show t, ϕ
    dpt = pressure_derivative(ϕ[1], Val(1), params.eos) #entropy
    dptt = pressure_derivative(ϕ[1], Val(2), params.eos)

    etaVisc = viscosity(ϕ[1], dpt, params.shear)
    tauShear = τ_shear(ϕ[1], dpt, params.shear)
    tauBulk = τ_bulk(ϕ[1], dpt, dptt, params.bulk)
    zetaBulkVisc = bulk_viscosity(ϕ[1], dpt, params.bulk)

    if free == true 
        thermo = hq_density(ϕ[1],ϕ[6];m=params.diffusion.mass)
            n = thermo.value
    dtn, dmn = thermo.gradient
    else 
        thermo = thermodynamic(ϕ[1],ϕ[6],params.eos.hadron_list)
            n = thermo.pressure
    dtn, dmn = thermo.pressure_derivative
    end
    

    dmn += 0.00001
    #dtn+=0.0001
    #dttn, dtdmn, dmmn = thermodynamic(ϕ[1],ϕ[6],HadronResonaceGasNew()).pressure_hessian.* fmGeV^3


    if free == true 

        Ds = diffusion(ϕ[1],n,params.diffusion)
        tauDiff= τ_diffusion(ϕ[1],params.diffusion)
    else 
        Ds = diffusion_hadron(ϕ[1],ϕ[6],params.eos,params.diffusion) #diffusion coefficient for hadrons
        tauDiff=τ_diffusion_hadron(ϕ[1],ϕ[6],params.eos,params.diffusion) #tau diffusion for hadrons
    end
    #@show t
    dmp = 0 #for now we don t have chemical potential in the eos
    dtdmp = 0
    dmdmp = 0

    #actually our equations don t depend on p: we can just put as entry dpt instead, in any case it will not be used (but in the future maybe it will be )

    #(At,Ax, source)=one_d_viscous_matrix(ϕ,t,X[1],dpt,dpt,dptt,dmp,dtdmp,dmdmp,zetaBulkVisc,etaVisc,tauShear,tauBulk,n,dtn,dmn,tauDiffusion,Ds)
    (At, Ax, source) = one_d_viscous_matrix_ruwen2(ϕ, t, X[1], dpt, dpt, dptt, zetaBulkVisc, etaVisc, tauShear, tauBulk, n, dtn, dmn, tauDiff, Ds)

    #println("Max(At)     = ", minimum(At),"     ",maximum(At))
    #println("Max(Ax)     = ", minimum(Ax),"     ",maximum(Ax))
    #println("Max(source) = ", minimum(source)," ",maximum(source))


    Ainv = inv(At)


    A_mul_B!(A_i[1], Ainv, Ax)


    jgemvavx!(Source, Ainv, source)



end


function one_d_viscous_matrix_ruwen2(X, tau, r, P_T, dP_dT, dP_dTdT, zetaBulkVisc, etaShearVisc, tauShear, tauBulk, n, dn_dT, dn_dmu, tauDiffusion, Ds)

    At = SMatrix{7,7}((X[2] .^2 .*dP_dT  .+ X[1] .*(1  .+ X[2] .^2) .*dP_dTdT,

X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ X[1] .*dP_dTdT),

0,

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

0,

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

0,

0,

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

tauDiffusion .*sqrt.(1  .+ X[2] .^2)))
    #########################################################################################################################################################################

    Ax = SMatrix{7,7}((X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ X[1] .*dP_dTdT),

(1  .+ X[2] .^2) .*dP_dT  .+ X[1] .*X[2] .^2 .*dP_dTdT,

0,

0,

0,

(3 .*tauShear .*tauBulk .*dP_dT .*(n .*X[1] .*X[2] .*(1  .+ X[2] .^2) .*dP_dTdT  .+ X[1] .*X[7] .*dP_dTdT  .- X[2] .*(1  .+ X[2] .^2) .*( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT) .*dn_dT)) ./((1  .+ X[2] .^2) .*(3 .*tauShear .*tauBulk .*X[1] .*X[2] .^2 .*dP_dT .^2  .+ X[1] .*(3 .*zetaBulkVisc .*tauShear .*X[2] .^2  .+ 4 .*etaShearVisc .*tauBulk .*X[2] .^2  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[3]  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[4]  .- 3 .*tauShear .*tauBulk .*X[5]  .- 3 .*tauShear .*tauBulk .*X[2] .^2 .*X[5]) .*dP_dTdT  .+ 3 .*tauShear .*tauBulk .*dP_dT .*( .-(X[2] .^2 .*X[3])  .- X[2] .^2 .*X[4]  .+ X[2] .^2 .*X[5]  .- X[1] .^2 .*dP_dTdT  .- X[1] .^2 .*X[2] .^2 .*dP_dTdT))),

(3 .*tauDiffusion .*tauShear .*tauBulk .*X[1] .*X[2] .*X[7] .*dP_dT .*dP_dTdT) ./( .-3 .*tauShear .*tauBulk .*X[1] .*X[2] .^2 .*dP_dT .^2  .- X[1] .*(3 .*zetaBulkVisc .*tauShear .*X[2] .^2  .+ 4 .*etaShearVisc .*tauBulk .*X[2] .^2  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[3]  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[4]  .- 3 .*tauShear .*tauBulk .*X[5]  .- 3 .*tauShear .*tauBulk .*X[2] .^2 .*X[5]) .*dP_dTdT  .+ 3 .*tauShear .*tauBulk .*dP_dT .*(X[2] .^2 .*X[3]  .+ X[2] .^2 .*X[4]  .- X[2] .^2 .*X[5]  .+ X[1] .^2 .*dP_dTdT  .+ X[1] .^2 .*X[2] .^2 .*dP_dTdT)),

((1  .+ 2 .*X[2] .^2) .*( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT)) ./sqrt.(1  .+ X[2] .^2),

2 .*X[2] .*( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT),

( .-2 .*etaShearVisc) ./(3. .*r .^2),

( .-2 .*etaShearVisc) ./(3. .*tau .^2),

zetaBulkVisc,

( .-(X[2] .*X[7] .*(3 .*tauShear .*tauBulk .*X[1] .*dP_dT .^2  .+ X[1] .*(3 .*zetaBulkVisc .*tauShear  .+ 4 .*etaShearVisc .*tauBulk  .+ 3 .*tauShear .*tauBulk .*X[3]  .+ 3 .*tauShear .*tauBulk .*X[4]  .- 3 .*tauShear .*tauBulk .*X[5]) .*dP_dTdT  .+ 3 .*tauShear .*tauBulk .*dP_dT .*( .-X[3]  .- X[4]  .+ X[5]  .- X[1] .^2 .*dP_dTdT)))  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT) .*( .-(n .*X[1] .*dP_dTdT)  .+ ( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT) .*dn_dT)) ./((1  .+ X[2] .^2) .*(3 .*tauShear .*tauBulk .*X[1] .*X[2] .^2 .*dP_dT .^2  .+ X[1] .*(3 .*zetaBulkVisc .*tauShear .*X[2] .^2  .+ 4 .*etaShearVisc .*tauBulk .*X[2] .^2  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[3]  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[4]  .- 3 .*tauShear .*tauBulk .*X[5]  .- 3 .*tauShear .*tauBulk .*X[2] .^2 .*X[5]) .*dP_dTdT  .+ 3 .*tauShear .*tauBulk .*dP_dT .*( .-(X[2] .^2 .*X[3])  .- X[2] .^2 .*X[4]  .+ X[2] .^2 .*X[5]  .- X[1] .^2 .*dP_dTdT  .- X[1] .^2 .*X[2] .^2 .*dP_dTdT))),

(tauDiffusion .*X[2] .^2 .*X[7] .*(3 .*tauShear .*tauBulk .*( .-X[3]  .- X[4]  .+ X[5]) .*dP_dT  .+ 3 .*tauShear .*tauBulk .*X[1] .*dP_dT .^2  .+ (3 .*zetaBulkVisc .*tauShear  .+ 4 .*etaShearVisc .*tauBulk) .*X[1] .*dP_dTdT)) ./((1  .+ X[2] .^2) .*(3 .*tauShear .*tauBulk .*X[1] .*X[2] .^2 .*dP_dT .^2  .+ X[1] .*(3 .*zetaBulkVisc .*tauShear .*X[2] .^2  .+ 4 .*etaShearVisc .*tauBulk .*X[2] .^2  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[3]  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[4]  .- 3 .*tauShear .*tauBulk .*X[5]  .- 3 .*tauShear .*tauBulk .*X[2] .^2 .*X[5]) .*dP_dTdT  .+ 3 .*tauShear .*tauBulk .*dP_dT .*( .-(X[2] .^2 .*X[3])  .- X[2] .^2 .*X[4]  .+ X[2] .^2 .*X[5]  .- X[1] .^2 .*dP_dTdT  .- X[1] .^2 .*X[2] .^2 .*dP_dTdT))),

 .-(X[2] .*sqrt.(1  .+ X[2] .^2)),

 .-1  .- X[2] .^2,

(tauShear .*X[2]) ./r .^2,

0,

0,

( .-3 .*tauShear .*tauBulk .*(n .*X[1] .*X[2] .*(1  .+ X[2] .^2) .*dP_dTdT  .+ X[1] .*X[7] .*dP_dTdT  .- X[2] .*(1  .+ X[2] .^2) .*( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT) .*dn_dT)) ./((1  .+ X[2] .^2) .*(3 .*tauShear .*tauBulk .*X[1] .*X[2] .^2 .*dP_dT .^2  .+ X[1] .*(3 .*zetaBulkVisc .*tauShear .*X[2] .^2  .+ 4 .*etaShearVisc .*tauBulk .*X[2] .^2  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[3]  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[4]  .- 3 .*tauShear .*tauBulk .*X[5]  .- 3 .*tauShear .*tauBulk .*X[2] .^2 .*X[5]) .*dP_dTdT  .+ 3 .*tauShear .*tauBulk .*dP_dT .*( .-(X[2] .^2 .*X[3])  .- X[2] .^2 .*X[4]  .+ X[2] .^2 .*X[5]  .- X[1] .^2 .*dP_dTdT  .- X[1] .^2 .*X[2] .^2 .*dP_dTdT))),

( .-3 .*tauDiffusion .*tauShear .*tauBulk .*X[1] .*X[2] .*X[7] .*dP_dTdT) ./( .-3 .*tauShear .*tauBulk .*X[1] .*X[2] .^2 .*dP_dT .^2  .- X[1] .*(3 .*zetaBulkVisc .*tauShear .*X[2] .^2  .+ 4 .*etaShearVisc .*tauBulk .*X[2] .^2  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[3]  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[4]  .- 3 .*tauShear .*tauBulk .*X[5]  .- 3 .*tauShear .*tauBulk .*X[2] .^2 .*X[5]) .*dP_dTdT  .+ 3 .*tauShear .*tauBulk .*dP_dT .*(X[2] .^2 .*X[3]  .+ X[2] .^2 .*X[4]  .- X[2] .^2 .*X[5]  .+ X[1] .^2 .*dP_dTdT  .+ X[1] .^2 .*X[2] .^2 .*dP_dTdT)),

 .-(X[2] .*sqrt.(1  .+ X[2] .^2)),

 .-1  .- X[2] .^2,

0,

(tauShear .*X[2]) ./tau .^2,

0,

( .-3 .*tauShear .*tauBulk .*(n .*X[1] .*X[2] .*(1  .+ X[2] .^2) .*dP_dTdT  .+ X[1] .*X[7] .*dP_dTdT  .- X[2] .*(1  .+ X[2] .^2) .*( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT) .*dn_dT)) ./((1  .+ X[2] .^2) .*(3 .*tauShear .*tauBulk .*X[1] .*X[2] .^2 .*dP_dT .^2  .+ X[1] .*(3 .*zetaBulkVisc .*tauShear .*X[2] .^2  .+ 4 .*etaShearVisc .*tauBulk .*X[2] .^2  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[3]  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[4]  .- 3 .*tauShear .*tauBulk .*X[5]  .- 3 .*tauShear .*tauBulk .*X[2] .^2 .*X[5]) .*dP_dTdT  .+ 3 .*tauShear .*tauBulk .*dP_dT .*( .-(X[2] .^2 .*X[3])  .- X[2] .^2 .*X[4]  .+ X[2] .^2 .*X[5]  .- X[1] .^2 .*dP_dTdT  .- X[1] .^2 .*X[2] .^2 .*dP_dTdT))),

( .-3 .*tauDiffusion .*tauShear .*tauBulk .*X[1] .*X[2] .*X[7] .*dP_dTdT) ./( .-3 .*tauShear .*tauBulk .*X[1] .*X[2] .^2 .*dP_dT .^2  .- X[1] .*(3 .*zetaBulkVisc .*tauShear .*X[2] .^2  .+ 4 .*etaShearVisc .*tauBulk .*X[2] .^2  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[3]  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[4]  .- 3 .*tauShear .*tauBulk .*X[5]  .- 3 .*tauShear .*tauBulk .*X[2] .^2 .*X[5]) .*dP_dTdT  .+ 3 .*tauShear .*tauBulk .*dP_dT .*(X[2] .^2 .*X[3]  .+ X[2] .^2 .*X[4]  .- X[2] .^2 .*X[5]  .+ X[1] .^2 .*dP_dTdT  .+ X[1] .^2 .*X[2] .^2 .*dP_dTdT)),

X[2] .*sqrt.(1  .+ X[2] .^2),

1  .+ X[2] .^2,

0,

0,

tauBulk .*X[2],

(3 .*tauShear .*tauBulk .*(n .*X[1] .*X[2] .*(1  .+ X[2] .^2) .*dP_dTdT  .+ X[1] .*X[7] .*dP_dTdT  .- X[2] .*(1  .+ X[2] .^2) .*( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT) .*dn_dT)) ./((1  .+ X[2] .^2) .*(3 .*tauShear .*tauBulk .*X[1] .*X[2] .^2 .*dP_dT .^2  .+ X[1] .*(3 .*zetaBulkVisc .*tauShear .*X[2] .^2  .+ 4 .*etaShearVisc .*tauBulk .*X[2] .^2  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[3]  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[4]  .- 3 .*tauShear .*tauBulk .*X[5]  .- 3 .*tauShear .*tauBulk .*X[2] .^2 .*X[5]) .*dP_dTdT  .+ 3 .*tauShear .*tauBulk .*dP_dT .*( .-(X[2] .^2 .*X[3])  .- X[2] .^2 .*X[4]  .+ X[2] .^2 .*X[5]  .- X[1] .^2 .*dP_dTdT  .- X[1] .^2 .*X[2] .^2 .*dP_dTdT))),

(3 .*tauDiffusion .*tauShear .*tauBulk .*X[1] .*X[2] .*X[7] .*dP_dTdT) ./( .-3 .*tauShear .*tauBulk .*X[1] .*X[2] .^2 .*dP_dT .^2  .- X[1] .*(3 .*zetaBulkVisc .*tauShear .*X[2] .^2  .+ 4 .*etaShearVisc .*tauBulk .*X[2] .^2  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[3]  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[4]  .- 3 .*tauShear .*tauBulk .*X[5]  .- 3 .*tauShear .*tauBulk .*X[2] .^2 .*X[5]) .*dP_dTdT  .+ 3 .*tauShear .*tauBulk .*dP_dT .*(X[2] .^2 .*X[3]  .+ X[2] .^2 .*X[4]  .- X[2] .^2 .*X[5]  .+ X[1] .^2 .*dP_dTdT  .+ X[1] .^2 .*X[2] .^2 .*dP_dTdT)),

0,

0,

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

tauDiffusion .*X[2]))
    #########################################################################################################################################################################

    source = SVector{7}(
        (( .-(X[2] .*X[3] .*(r .*X[2]  .+ sqrt.(1  .+ X[2] .^2) .*tau))  .+ X[4] .*(r  .- r .*X[2] .^2  .- X[2] .*sqrt.(1  .+ X[2] .^2) .*tau)  .+ (r  .+ r .*X[2] .^2  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau) .*(X[5]  .+ X[1] .*dP_dT)) ./(r .*tau),

( .-(X[4] .*(r .*X[2] .*sqrt.(1  .+ X[2] .^2)  .+ tau  .+ X[2] .^2 .*tau))  .- X[3] .*(r .*X[2] .*sqrt.(1  .+ X[2] .^2)  .+ 2 .*tau  .+ X[2] .^2 .*tau)  .+ X[2] .*(r .*sqrt.(1  .+ X[2] .^2)  .+ X[2] .*tau) .*(X[5]  .+ X[1] .*dP_dT)) ./(r .*tau),

( .-2 .*etaShearVisc .*r .*sqrt.(1  .+ X[2] .^2)  .+ 4 .*etaShearVisc .*X[2] .*tau  .+ 3 .*r .*X[3] .*tau) ./(3. .*r .^3 .*tau),

(4 .*etaShearVisc .*r .*sqrt.(1  .+ X[2] .^2)  .- 2 .*etaShearVisc .*X[2] .*tau  .+ 3 .*r .*X[4] .*tau) ./(3. .*r .*tau .^3),

(zetaBulkVisc .*X[2]) ./r  .+ X[5]  .+ (zetaBulkVisc .*sqrt.(1  .+ X[2] .^2)) ./tau,

(X[7] .*(3 .*tauShear .*tauBulk .*X[1] .*X[2] .*(r .*( .-1  .+ X[2] .^4)  .+ X[2] .^3 .*sqrt.(1  .+ X[2] .^2) .*tau) .*dP_dT .^2  .+ X[1] .*(r .*X[2] .*(1  .+ X[2] .^2) .*(3 .*zetaBulkVisc .*tauShear .*( .-1  .+ X[2] .^2)  .+ etaShearVisc .*tauBulk .*(2  .+ 4 .*X[2] .^2)  .- 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*X[5])  .+ sqrt.(1  .+ X[2] .^2) .*( .-3 .*tauShear .*tauBulk .*X[5]  .- 3 .*tauShear .*r .*X[2] .*X[5]  .+ 6 .*tauBulk .*X[2] .^2 .*(etaShearVisc  .- tauShear .*X[5])  .+ X[2] .^4 .*(3 .*zetaBulkVisc .*tauShear  .+ 4 .*etaShearVisc .*tauBulk  .- 3 .*tauShear .*tauBulk .*X[5])) .*tau  .+ 3 .*tauBulk .*X[3] .*(tauShear .*r .*X[2] .^3 .*(1  .+ X[2] .^2)  .- tauShear .*sqrt.(1  .+ X[2] .^2) .*tau  .+ r .*X[2] .*sqrt.(1  .+ X[2] .^2) .*tau  .+ tauShear .*X[2] .^4 .*sqrt.(1  .+ X[2] .^2) .*tau)  .+ 3 .*tauBulk .*X[2] .*X[4] .*(r .*sqrt.(1  .+ X[2] .^2) .*tau  .+ tauShear .*(1  .+ X[2] .^2) .*(r .*( .-1  .+ X[2] .^2)  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau))) .*dP_dTdT  .+ 3 .*tauShear .*tauBulk .*dP_dT .*(r .*X[2] .*( .-1  .+ X[2] .^4) .*X[5]  .+ X[2] .^4 .*sqrt.(1  .+ X[2] .^2) .*X[5] .*tau  .- X[2] .*(1  .+ X[2] .^2) .*X[4] .*(r  .+ r .*X[2] .^2  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau)  .- X[2] .^2 .*X[3] .*(2 .*sqrt.(1  .+ X[2] .^2) .*tau  .+ X[2] .*(r  .+ r .*X[2] .^2  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau))  .- X[1] .^2 .*(1  .+ X[2] .^2) .^2 .*(r .*X[2]  .+ sqrt.(1  .+ X[2] .^2) .*tau) .*dP_dTdT))  .- 3 .*n .*(1  .+ X[2] .^2) .*(X[1] .*(tauShear .*tauBulk .*r .*X[5]  .+ tauBulk .*r .*X[2] .^4 .*( .-2 .*etaShearVisc  .+ tauShear .*X[5])  .+ tauShear .*tauBulk .*X[2] .*sqrt.(1  .+ X[2] .^2) .*X[5] .*tau  .+ tauBulk .*X[2] .^3 .*sqrt.(1  .+ X[2] .^2) .*( .-2 .*etaShearVisc  .+ tauShear .*X[5]) .*tau  .+ r .*X[2] .^2 .*( .-2 .*etaShearVisc .*tauBulk  .+ tauShear .*X[5] .*(2 .*tauBulk  .+ sqrt.(1  .+ X[2] .^2) .*tau))  .+ tauShear .*tauBulk .*X[1] .*(1  .+ X[2] .^2) .*(r  .+ r .*X[2] .^2  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau) .*dP_dT) .*dP_dTdT  .+ tauBulk .*X[4] .*(tauShear .*X[2] .^2 .*(2 .*r .*(1  .+ X[2] .^2)  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau) .*dP_dT  .+ r .*X[1] .*(tauShear .*( .-1  .+ X[2] .^4)  .- X[2] .^2 .*sqrt.(1  .+ X[2] .^2) .*tau) .*dP_dTdT)  .+ tauBulk .*X[3] .*(tauShear .*X[2] .^2 .*(r  .+ r .*X[2] .^2  .+ 2 .*X[2] .*sqrt.(1  .+ X[2] .^2) .*tau) .*dP_dT  .- X[1] .*(r .*X[2] .^2 .*sqrt.(1  .+ X[2] .^2) .*tau  .+ tauShear .*(1  .+ X[2] .^2) .*(r  .- X[2] .*sqrt.(1  .+ X[2] .^2) .*tau)) .*dP_dTdT))  .- (1  .+ X[2] .^2) .*(3 .*tauBulk .*X[2] .*X[3] .^2 .*(tauShear .*r .*X[2] .^3  .+ 3 .*tauShear .*sqrt.(1  .+ X[2] .^2) .*tau  .+ 3 .*tauShear .*X[2] .^2 .*sqrt.(1  .+ X[2] .^2) .*tau  .+ r .*X[2] .*(tauShear  .- sqrt.(1  .+ X[2] .^2) .*tau))  .+ 3 .*tauBulk .*X[4] .^2 .*( .-(r .*X[2] .^2 .*sqrt.(1  .+ X[2] .^2) .*tau)  .+ tauShear .*(1  .+ X[2] .^2) .*(r  .+ 3 .*r .*X[2] .^2  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau))  .- 3 .*(X[5]  .+ X[1] .*dP_dT) .*(tauShear .*tauBulk .*r .*X[5]  .+ tauBulk .*r .*X[2] .^4 .*( .-2 .*etaShearVisc  .+ tauShear .*X[5])  .+ tauShear .*tauBulk .*X[2] .*sqrt.(1  .+ X[2] .^2) .*X[5] .*tau  .+ tauBulk .*X[2] .^3 .*sqrt.(1  .+ X[2] .^2) .*( .-2 .*etaShearVisc  .+ tauShear .*X[5]) .*tau  .+ r .*X[2] .^2 .*( .-2 .*etaShearVisc .*tauBulk  .+ tauShear .*X[5] .*(2 .*tauBulk  .+ sqrt.(1  .+ X[2] .^2) .*tau))  .+ tauShear .*tauBulk .*X[1] .*(1  .+ X[2] .^2) .*(r  .+ r .*X[2] .^2  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau) .*dP_dT)  .+ X[3] .*(r .*(1  .+ X[2] .^2) .*((3 .*zetaBulkVisc .*tauShear  .- 2 .*etaShearVisc .*tauBulk) .*X[2] .^2  .+ 3 .*tauShear .*tauBulk .*X[5])  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*( .-6 .*tauShear .*tauBulk .*X[5]  .+ 3 .*(tauShear  .+ tauBulk) .*r .*X[2] .*X[5]  .+ 2 .*X[2] .^2 .*(3 .*zetaBulkVisc .*tauShear  .+ etaShearVisc .*tauBulk  .- 3 .*tauShear .*tauBulk .*X[5])) .*tau  .+ 3 .*tauBulk .*X[4] .*( .-2 .*r .*X[2] .^2 .*sqrt.(1  .+ X[2] .^2) .*tau  .+ tauShear .*(1  .+ X[2] .^2) .*(r  .+ 4 .*r .*X[2] .^2  .+ 4 .*X[2] .*sqrt.(1  .+ X[2] .^2) .*tau))  .+ 3 .*tauBulk .*X[1] .*(r .*X[2] .^2 .*sqrt.(1  .+ X[2] .^2) .*tau  .+ tauShear .*(1  .+ X[2] .^2) .*(r  .- 2 .*X[2] .*sqrt.(1  .+ X[2] .^2) .*tau)) .*dP_dT)  .+ X[2] .^2 .*X[4] .*(3 .*zetaBulkVisc .*tauShear .*X[2] .*sqrt.(1  .+ X[2] .^2) .*tau  .- 2 .*etaShearVisc .*tauBulk .*X[2] .*sqrt.(1  .+ X[2] .^2) .*tau  .+ 3 .*tauShear .*r .*sqrt.(1  .+ X[2] .^2) .*X[5] .*tau  .+ 3 .*tauBulk .*r .*sqrt.(1  .+ X[2] .^2) .*X[5] .*tau  .+ 3 .*tauBulk .*r .*X[1] .*sqrt.(1  .+ X[2] .^2) .*tau .*dP_dT  .+ 2 .*r .*(1  .+ X[2] .^2) .*(3 .*zetaBulkVisc .*tauShear  .+ etaShearVisc .*tauBulk  .- 3 .*tauShear .*tauBulk .*X[5]  .- 3 .*tauShear .*tauBulk .*X[1] .*dP_dT))) .*dn_dT) ./(r .*(1  .+ X[2] .^2) .^1.5 .*tau .*(3 .*tauShear .*tauBulk .*X[2] .^2 .*dP_dT .*( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT)  .+ X[1] .*( .-3 .*tauShear .*tauBulk .*X[5]  .+ X[2] .^2 .*(3 .*zetaBulkVisc .*tauShear  .+ 4 .*etaShearVisc .*tauBulk  .- 3 .*tauShear .*tauBulk .*X[5])  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*(X[3]  .+ X[4]  .- X[1] .*dP_dT)) .*dP_dTdT)),

 .-((X[7] .*( .-3 .*tauShear .*tauBulk .*X[1] .*X[2] .^2 .*(r .*sqrt.(1  .+ X[2] .^2) .*tau  .+ tauDiffusion .*(r  .+ r .*X[2] .^2  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau)) .*dP_dT .^2  .- X[1] .*(tauDiffusion .*(3 .*zetaBulkVisc .*tauShear  .- 2 .*etaShearVisc .*tauBulk) .*r .*X[2] .^2 .*(1  .+ X[2] .^2)  .+ sqrt.(1  .+ X[2] .^2) .*(tauDiffusion .*(3 .*zetaBulkVisc .*tauShear  .- 2 .*etaShearVisc .*tauBulk) .*X[2] .^3  .- 3 .*tauShear .*tauBulk .*r .*X[5]  .+ r .*X[2] .^2 .*(4 .*etaShearVisc .*tauBulk  .+ 3 .*tauShear .*(zetaBulkVisc  .+ tauDiffusion .*X[5]  .- tauBulk .*X[5]))) .*tau  .+ 3 .*tauBulk .*X[3] .*(tauDiffusion .*tauShear .*r .*X[2] .^4  .+ tauShear .*r .*sqrt.(1  .+ X[2] .^2) .*tau  .+ 2 .*tauDiffusion .*tauShear .*X[2] .*sqrt.(1  .+ X[2] .^2) .*tau  .+ 2 .*tauDiffusion .*tauShear .*X[2] .^3 .*sqrt.(1  .+ X[2] .^2) .*tau  .+ r .*X[2] .^2 .*(tauDiffusion .*tauShear  .- tauDiffusion .*sqrt.(1  .+ X[2] .^2) .*tau  .+ tauShear .*sqrt.(1  .+ X[2] .^2) .*tau))  .+ 3 .*tauBulk .*X[4] .*(2 .*tauDiffusion .*tauShear .*r .*X[2] .^4  .+ tauShear .*r .*sqrt.(1  .+ X[2] .^2) .*tau  .+ tauDiffusion .*tauShear .*X[2] .*sqrt.(1  .+ X[2] .^2) .*tau  .+ tauDiffusion .*tauShear .*X[2] .^3 .*sqrt.(1  .+ X[2] .^2) .*tau  .+ r .*X[2] .^2 .*(2 .*tauDiffusion .*tauShear  .- tauDiffusion .*sqrt.(1  .+ X[2] .^2) .*tau  .+ tauShear .*sqrt.(1  .+ X[2] .^2) .*tau))) .*dP_dTdT  .- 3 .*tauShear .*tauBulk .*dP_dT .*(tauDiffusion .*r .*X[2] .^2 .*(1  .+ X[2] .^2) .*X[5]  .+ r .*X[2] .^2 .*X[4] .*(tauDiffusion  .+ tauDiffusion .*X[2] .^2  .- sqrt.(1  .+ X[2] .^2) .*tau)  .+ sqrt.(1  .+ X[2] .^2) .*tau .*(X[2] .^2 .*(( .-r  .+ tauDiffusion .*X[2]) .*X[3]  .+ (r  .+ tauDiffusion .*X[2]) .*X[5])  .- r .*X[1] .^2 .*(1  .+ X[2] .^2) .*dP_dTdT)))) ./(r .*sqrt.(1  .+ X[2] .^2) .*tau .*(3 .*tauShear .*tauBulk .*X[2] .^2 .*dP_dT .*( .-X[3]  .- X[4]  .+ X[5]  .+ X[1] .*dP_dT)  .+ X[1] .*( .-3 .*tauShear .*tauBulk .*X[5]  .+ X[2] .^2 .*(3 .*zetaBulkVisc .*tauShear  .+ 4 .*etaShearVisc .*tauBulk  .- 3 .*tauShear .*tauBulk .*X[5])  .+ 3 .*tauShear .*tauBulk .*(1  .+ X[2] .^2) .*(X[3]  .+ X[4]  .- X[1] .*dP_dT)) .*dP_dTdT))))
    )
         return (At,Ax, source)
end

