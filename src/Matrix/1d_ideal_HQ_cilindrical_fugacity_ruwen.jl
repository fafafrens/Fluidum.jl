# the convention here are T, ur,  \[Pi]phiphi, \[Pi]etaeta, \[Pi]B, α, nur

function matrix1d_ideal_HQ_ruwen!(A_i, Source, ϕ, t, X, params;free=true)

    #@show t, ϕ
    dpt = pressure_derivative(ϕ[1], Val(1), params.eos) #entropy
    dptt = pressure_derivative(ϕ[1], Val(2), params.eos)

    if free == true 
        thermo = hq_density(ϕ[1],ϕ[3];m=params.diffusion.mass)
            n = thermo.value
    dtn, dmn = thermo.gradient
    else 
        thermo = thermodynamic(ϕ[1],ϕ[3],params.eos.hadron_list)
            n = thermo.pressure
    dtn, dmn = thermo.pressure_derivative
    end

    dmn += 1e-6
    dtn += 1e-6
    #actually our equations don t depend on p: we can just put as entry dpt instead, in any case it will not be used (but in the future maybe it will be )

    #(At,Ax, source)=one_d_viscous_matrix(ϕ,t,X[1],dpt,dpt,dptt,dmp,dtdmp,dmdmp,zetaBulkVisc,etaVisc,tauShear,tauBulk,n,dtn,dmn,tauDiffusion,Ds)
    (At, Ax, source) = one_d_ideal_matrix_ruwen(ϕ, t, X[1], dpt, dpt, dptt, n, dtn, dmn)

    #println("Max(At)     = ", minimum(At),"     ",maximum(At))
    #println("Max(Ax)     = ", minimum(Ax),"     ",maximum(Ax))
    #println("Max(source) = ", minimum(source)," ",maximum(source))


    Ainv = inv(At)


    A_mul_B!(A_i[1], Ainv, Ax)


    jgemvavx!(Source, Ainv, source)



end


function one_d_ideal_matrix_ruwen(X, tau, r, P_T, dP_dT, dP_dTdT, n, dn_dT, dn_dmu)

    At = SMatrix{3,3}(X[2] .^2 .*dP_dT  .+ X[1] .*(1  .+ X[2] .^2) .*dP_dTdT,

X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ X[1] .*dP_dTdT),

0,

2 .*X[1] .*X[2] .*dP_dT,

(X[1] .*(1  .+ 2 .*X[2] .^2) .*dP_dT) ./sqrt.(1  .+ X[2] .^2),

0,

0,

0,

sqrt.(1  .+ X[2] .^2) .*dn_dmu)
    #########################################################################################################################################################################

    Ax = SMatrix{3,3}(X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ X[1] .*dP_dTdT),

(1  .+ X[2] .^2) .*dP_dT  .+ X[1] .*X[2] .^2 .*dP_dTdT,

(X[2] .*(n .*dP_dTdT  .- dP_dT .*dn_dT)) ./(X[2] .^2 .*dP_dT  .- X[1] .*(1  .+ X[2] .^2) .*dP_dTdT),

(X[1] .*(1  .+ 2 .*X[2] .^2) .*dP_dT) ./sqrt.(1  .+ X[2] .^2),

2 .*X[1] .*X[2] .*dP_dT,

( .-(n .*X[1] .*dP_dTdT)  .+ X[1] .*dP_dT .*dn_dT) ./(X[2] .^2 .*dP_dT  .- X[1] .*(1  .+ X[2] .^2) .*dP_dTdT),

0,

0,

X[2] .*dn_dmu)
    #########################################################################################################################################################################

    source = SVector{3}((X[1] .*sqrt.(1  .+ X[2] .^2) .*(r .*sqrt.(1  .+ X[2] .^2)  .+ X[2] .*tau) .*dP_dT) ./(r .*tau),

(X[1] .*X[2] .*(r  .+ r .*X[2] .^2  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau) .*dP_dT) ./(r .*sqrt.(1  .+ X[2] .^2) .*tau),

(X[1] .*(1  .+ X[2] .^2) .*(r .*sqrt.(1  .+ X[2] .^2)  .+ X[2] .*tau) .*( .-(n .*dP_dTdT)  .+ dP_dT .*dn_dT)) ./(r .*tau .*(X[2] .^2 .*dP_dT  .- X[1] .*(1  .+ X[2] .^2) .*dP_dTdT)))
         return (At,Ax, source)
end

