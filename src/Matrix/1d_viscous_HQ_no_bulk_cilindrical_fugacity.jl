# the convention here are T, ur,  \[Pi]phiphi, \[Pi]etaeta, \[Pi]B, α, nur

function matrix1d_visc_no_bulk_HQ!(A_i, Source, ϕ, t, X, params;free=true)

    #@show t, ϕ
    dpt = pressure_derivative(ϕ[1], Val(1), params.eos) #entropy
    dptt = pressure_derivative(ϕ[1], Val(2), params.eos)

    etaVisc = viscosity(ϕ[1], dpt, params.shear)
    tauShear = τ_shear(ϕ[1], dpt, params.shear)
    #tauBulk = τ_bulk(ϕ[1], dpt, dptt, params.bulk)
    #zetaBulkVisc = bulk_viscosity(ϕ[1], dpt, params.bulk)

    if free == true 
        thermo = hq_density(ϕ[1],ϕ[5];m=params.diffusion.mass)
            n = thermo.value
            dtn, dmn = thermo.gradient
    else 
        thermo = thermodynamic(ϕ[1],ϕ[5],params.eos.hadron_list)
        n = thermo.pressure
        dtn, dmn = thermo.pressure_derivative
    end
    

    dmn += 0.00001
    dtn += 0.00001
    #dttn, dtdmn, dmmn = thermodynamic(ϕ[1],ϕ[6],HadronResonaceGasNew()).pressure_hessian.* fmGeV^3


    if free == true 

        Ds = diffusion(ϕ[1],n,params.diffusion)
        tauDiff= τ_diffusion(ϕ[1],params.diffusion)
    else 
        Ds = diffusion_hadron(ϕ[1],ϕ[5],params.eos,params.diffusion) #diffusion coefficient for hadrons
        tauDiff=τ_diffusion_hadron(ϕ[1],ϕ[5],params.eos,params.diffusion) #tau diffusion for hadrons
    end


    #@debug "HQ diag" diag_HQ

    #@show t
    dmp = 0 #for now we don t have chemical potential in the eos
    dtdmp = 0
    dmdmp = 0

    #actually our equations don t depend on p: we can just put as entry dpt instead, in any case it will not be used (but in the future maybe it will be )

    #(At,Ax, source)=one_d_viscous_matrix(ϕ,t,X[1],dpt,dpt,dptt,dmp,dtdmp,dmdmp,zetaBulkVisc,etaVisc,tauShear,tauBulk,n,dtn,dmn,tauDiffusion,Ds)
    (At, Ax, source) = one_d_viscous_matrix9(ϕ, t, X[1], dpt, dpt, dptt, dpt, etaVisc, tauShear, dpt, n, dtn, dmn, tauDiff, Ds)

    #println("Max(At)     = ", minimum(At),"     ",maximum(At))
    #println("Max(Ax)     = ", minimum(Ax),"     ",maximum(Ax))
    #println("Max(source) = ", minimum(source)," ",maximum(source))


    Ainv = inv(At)


    A_mul_B!(A_i[1], Ainv, Ax)


    jgemvavx!(Source, Ainv, source)



end


function one_d_viscous_matrix9(X, tau, r, P_T, dP_dT, dP_dTdT, zetaBulkVisc, etaShearVisc, tauShear, tauBulk, n, dn_dT, dn_dmu, tauDiffusion, Ds)

    At = SMatrix{6,6}(
(X[2] .^2 .*dP_dT  .+ X[1] .*(1  .+ X[2] .^2) .*dP_dTdT,

X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ X[1] .*dP_dTdT),

0,

0,

sqrt.(1  .+ X[2] .^2) .*dn_dT,

0,

 .-2 .*X[2] .*(X[3]  .+ X[4]  .- X[1] .*dP_dT),

 .-(((1  .+ 2 .*X[2] .^2) .*(X[3]  .+ X[4]  .- X[1] .*dP_dT)) ./sqrt.(1  .+ X[2] .^2)),

( .-2 .*etaShearVisc .*X[2]) ./(3. .*r .^2 .*sqrt.(1  .+ X[2] .^2)),

( .-2 .*etaShearVisc .*X[2]) ./(3. .*sqrt.(1  .+ X[2] .^2) .*tau .^2),

(n .*(X[2]  .+ X[2] .^3)  .+ X[6]) ./(1  .+ X[2] .^2) .^1.5,

 .-((tauDiffusion .*X[2] .*X[6]) ./sqrt.(1  .+ X[2] .^2)),

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

X[2] ./sqrt.(1  .+ X[2] .^2),

tauDiffusion .*sqrt.(1  .+ X[2] .^2))
)
    #########################################################################################################################################################################

    Ax = SMatrix{6,6}(
(X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ X[1] .*dP_dTdT),

(1  .+ X[2] .^2) .*dP_dT  .+ X[1] .*X[2] .^2 .*dP_dTdT,

0,

0,

X[2] .*dn_dT,

0,

 .-(((1  .+ 2 .*X[2] .^2) .*(X[3]  .+ X[4]  .- X[1] .*dP_dT)) ./sqrt.(1  .+ X[2] .^2)),

 .-2 .*X[2] .*(X[3]  .+ X[4]  .- X[1] .*dP_dT),

( .-2 .*etaShearVisc) ./(3. .*r .^2),

( .-2 .*etaShearVisc) ./(3. .*tau .^2),

n,

 .-((tauDiffusion .*X[2] .^2 .*X[6]) ./(1  .+ X[2] .^2)),

 .-(X[2] .*sqrt.(1  .+ X[2] .^2)),

 .-1  .- X[2] .^2,

(tauShear .*X[2]) ./r .^2,

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

0,

0,

0,

X[2] .*dn_dmu,

Ds .*n .*(1  .+ X[2] .^2),

0,

0,

0,

0,

1,

tauDiffusion .*X[2])
    )
    #########################################################################################################################################################################

    source = SVector{6}(
(( .-(X[2] .*X[3] .*(r .*X[2]  .+ sqrt.(1  .+ X[2] .^2) .*tau))  .+ X[4] .*(r  .- r .*X[2] .^2  .- X[2] .*sqrt.(1  .+ X[2] .^2) .*tau)  .+ X[1] .*(r  .+ r .*X[2] .^2  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau) .*dP_dT) ./(r .*tau),

( .-(X[4] .*(r .*X[2] .*sqrt.(1  .+ X[2] .^2)  .+ tau  .+ X[2] .^2 .*tau))  .- X[3] .*(r .*X[2] .*sqrt.(1  .+ X[2] .^2)  .+ 2 .*tau  .+ X[2] .^2 .*tau)  .+ X[1] .*X[2] .*(r .*sqrt.(1  .+ X[2] .^2)  .+ X[2] .*tau) .*dP_dT) ./(r .*tau),

( .-2 .*etaShearVisc .*r .*sqrt.(1  .+ X[2] .^2)  .+ 4 .*etaShearVisc .*X[2] .*tau  .+ 3 .*r .*X[3] .*tau) ./(3. .*r .^3 .*tau),

(4 .*etaShearVisc .*r .*sqrt.(1  .+ X[2] .^2)  .- 2 .*etaShearVisc .*X[2] .*tau  .+ 3 .*r .*X[4] .*tau) ./(3. .*r .*tau .^3),

X[6] .*(1 ./r  .+ X[2] ./(sqrt.(1  .+ X[2] .^2) .*tau))  .+ n .*(X[2] ./r  .+ sqrt.(1  .+ X[2] .^2) ./tau),

X[6])
    )
         return (At,Ax, source)
end

