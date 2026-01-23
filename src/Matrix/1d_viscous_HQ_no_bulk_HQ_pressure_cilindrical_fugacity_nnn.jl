# the convention here are T, ur,  \[Pi]phiphi, \[Pi]etaeta, \[Pi]B, α, nur

function matrix1d_visc_no_bulk_pressure_HQ!(A_i, Source, ϕ, t, X, params;free=true)

    T = ϕ[1]
    α = ϕ[5]

    # --- 1. Check basic validity -----------------------------------
    if !isfinite(T) || !isfinite(α)
        @error "Non-finite state detected" T=T α=α t=t X=X
    end

    if T <= 0
        @error "NEGATIVE TEMPERATURE detected" T=T α=α t=t X=X
        T = max(T, 1e-9)   # keep code running
    end

    #if abs(α) > 50
    #    @warn "Fugacity α extremely large. Expect overflow in exp(α)" α=α t=t
    #end

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
    dn_dTdT, dn_dTdmu,_ = thermo.pressure_hessian
    #dttn, dtdmn, dmmn = thermodynamic(ϕ[1],ϕ[6],HadronResonaceGasNew()).pressure_hessian.* fmGeV^3
    #dn_dTdmu = 0
    if free == true 
        Ds = diffusion(ϕ[1],n,params.diffusion)
        tauDiff= τ_diffusion(ϕ[1],params.diffusion)
    else 
        Ds = diffusion_hadron(ϕ[1],ϕ[5],params.eos,params.diffusion) #diffusion coefficient for hadrons
        tauDiff=τ_diffusion_hadron(ϕ[1],ϕ[5],params.eos,params.diffusion) #tau diffusion for hadrons
    end

    pressure_field = hq_pressure(ϕ[1],ϕ[5];m=params.diffusion.mass)
    Phq_Tmu = pressure_field.pressure
    dPhq_dT, dPhq_dmu = pressure_field.pressure_derivative
    #dPhq_dTdT,dPhq_dTdmu,ddp_dmudmu_HQ = pressure_field.pressure_hessian
    #dPhq_dT = 0.
    #dPhq_dmu = 0.
    #dPhq_dTdT = 0.
    #dPhq_dTdmu = 0.

    (At, Ax, source) = one_d_viscous_matrix10(ϕ, t, X[1], nothing, dpt, dptt,Phq_Tmu,dPhq_dT,dPhq_dmu,nothing, nothing, nothing, etaVisc, tauShear, nothing, n, dtn,dmn,dn_dTdT,dn_dTdmu, tauDiff, Ds)

    #debug_negative_temperature(ϕ, X, dPhq_dmu, dPhq_dTdmu)
    #println("Max(At)     = ", minimum(At),"     ",maximum(At))
    #println("Max(Ax)     = ", minimum(Ax),"     ",maximum(Ax))
    #println("Max(source) = ", minimum(source)," ",maximum(source))


    Ainv = inv(At)


    A_mul_B!(A_i[1], Ainv, Ax)


    jgemvavx!(Source, Ainv, source)



end


function one_d_viscous_matrix10(X, tau, r, P_T, dP_dT, dP_dTdT,Phq_Tmu,dPhq_dT, dPhq_dmu,dPhq_dTdT, dPhq_dTdmu, zetaBulkVisc, etaShearVisc, tauShear, tauBulk, n, dn_dT,dn_dmu,dn_dTdT,dn_dTdmu, tauDiffusion, Ds)

    At = SMatrix{6,6}(
    (X[2] .^2 .*dP_dT  .+ X[2] .^2 .*dPhq_dT  .+ X[1] .*(1  .+ X[2] .^2) .*(dP_dTdT  .+ 2 .*dn_dT  .+ X[1] .*dn_dTdT),

X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ dPhq_dT  .+ X[1] .*(dP_dTdT  .+ 2 .*dn_dT  .+ X[1] .*dn_dTdT)),

0,

0,

sqrt.(1  .+ X[2] .^2) .*dn_dT,

0,

2 .*X[2] .*(Phq_Tmu  .- X[3]  .- X[4]  .+ X[1] .*dP_dT  .+ X[1] .^2 .*dn_dT),

((1  .+ 2 .*X[2] .^2) .*(Phq_Tmu  .- X[3]  .- X[4]  .+ X[1] .*dP_dT  .+ X[1] .^2 .*dn_dT)) ./sqrt.(1  .+ X[2] .^2),

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

X[2] .^2 .*dPhq_dmu  .+ X[1] .^2 .*(1  .+ X[2] .^2) .*dn_dTdmu,

X[2] .*sqrt.(1  .+ X[2] .^2) .*(dPhq_dmu  .+ X[1] .^2 .*dn_dTdmu),

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
(X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ dPhq_dT  .+ X[1] .*(dP_dTdT  .+ 2 .*dn_dT  .+ X[1] .*dn_dTdT)),

(1  .+ X[2] .^2) .*dP_dT  .+ (1  .+ X[2] .^2) .*dPhq_dT  .+ X[1] .*X[2] .^2 .*(dP_dTdT  .+ 2 .*dn_dT  .+ X[1] .*dn_dTdT),

0,

0,

X[2] .*dn_dT,

0,

((1  .+ 2 .*X[2] .^2) .*(Phq_Tmu  .- X[3]  .- X[4]  .+ X[1] .*dP_dT  .+ X[1] .^2 .*dn_dT)) ./sqrt.(1  .+ X[2] .^2),

2 .*X[2] .*(Phq_Tmu  .- X[3]  .- X[4]  .+ X[1] .*dP_dT  .+ X[1] .^2 .*dn_dT),

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

X[2] .*sqrt.(1  .+ X[2] .^2) .*(dPhq_dmu  .+ X[1] .^2 .*dn_dTdmu),

(1  .+ X[2] .^2) .*dPhq_dmu  .+ X[1] .^2 .*X[2] .^2 .*dn_dTdmu,

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
((r .*X[4]  .- r .*X[2] .^2 .*X[4]  .- X[2] .*sqrt.(1  .+ X[2] .^2) .*X[4] .*tau  .- X[2] .*X[3] .*(r .*X[2]  .+ sqrt.(1  .+ X[2] .^2) .*tau)  .+ Phq_Tmu .*(r  .+ r .*X[2] .^2  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau)  .+ r .*X[1] .*dP_dT  .+ r .*X[1] .*X[2] .^2 .*dP_dT  .+ X[1] .*X[2] .*sqrt.(1  .+ X[2] .^2) .*tau .*dP_dT  .+ r .*X[1] .^2 .*dn_dT  .+ r .*X[1] .^2 .*X[2] .^2 .*dn_dT  .+ X[1] .^2 .*X[2] .*sqrt.(1  .+ X[2] .^2) .*tau .*dn_dT) ./(r .*tau),

( .-(r .*X[2] .*sqrt.(1  .+ X[2] .^2) .*X[4])  .- X[4] .*tau  .+ Phq_Tmu .*X[2] .*(r .*sqrt.(1  .+ X[2] .^2)  .+ X[2] .*tau)  .- X[3] .*(r .*X[2] .*sqrt.(1  .+ X[2] .^2)  .+ 2 .*tau  .+ X[2] .^2 .*tau)  .+ r .*X[1] .*X[2] .*sqrt.(1  .+ X[2] .^2) .*dP_dT  .+ r .*X[1] .^2 .*X[2] .*sqrt.(1  .+ X[2] .^2) .*dn_dT  .+ X[2] .^2 .*tau .*( .-X[4]  .+ X[1] .*(dP_dT  .+ X[1] .*dn_dT))) ./(r .*tau),

( .-2 .*etaShearVisc .*r .*sqrt.(1  .+ X[2] .^2)  .+ 4 .*etaShearVisc .*X[2] .*tau  .+ 3 .*r .*X[3] .*tau) ./(3. .*r .^3 .*tau),

(4 .*etaShearVisc .*r .*sqrt.(1  .+ X[2] .^2)  .- 2 .*etaShearVisc .*X[2] .*tau  .+ 3 .*r .*X[4] .*tau) ./(3. .*r .*tau .^3),

X[6] .*(1 ./r  .+ X[2] ./(sqrt.(1  .+ X[2] .^2) .*tau))  .+ n .*(X[2] ./r  .+ sqrt.(1  .+ X[2] .^2) ./tau),

X[6])
    )
         return (At,Ax, source)
end

