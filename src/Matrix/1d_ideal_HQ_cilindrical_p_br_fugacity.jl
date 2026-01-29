#
# The matrices are based on the notebook: ideal_r_tau_br_p_fugacity.nb
# The fields are stored in the order : T, ur, α
# HQ back reaction and pressure feedback included
#

function matrix1d_ideal_HQ_br_p_fugacity!(A_i, Source, ϕ, t, X, params;free=true)

    #@show t, ϕ
    dP_dT = pressure_derivative(ϕ[1], Val(1), params.eos) #entropy
    dP_dTdT = pressure_derivative(ϕ[1], Val(2), params.eos)



    if free == true 
        therm = hq_pressure(ϕ[1], ϕ[3]; m=params.diffusion.mass)    
        dPhq_dT, dPhq_dalpha = therm.pressure_derivative
        ddPhq_dTdT, ddPhq_dTdalpha, ddPhq_dalphadalpha = therm.pressure_hessian

        thermo = hq_density(ϕ[1], ϕ[3]; m=params.diffusion.mass)
        n = thermo.value
        dn_dT, dn_dmu = thermo.gradient
    else
       @warn "not free ideal not implemented."
    end

    dmn_eps = let s = get(ENV, "FLUIDUM_DMN_EPS", "1e-4")
    v = tryparse(Float64, s)
    v === nothing ? 1e-4 : v
    end
    dtn_eps = let s = get(ENV, "FLUIDUM_DTN_EPS", "1e-4")
        v = tryparse(Float64, s)
        v === nothing ? 1e-4 : v
    end

    #dn_dmu += dmn_eps
    #dn_dT += dtn_eps

    (At, Ax, source) = one_d_ideal_matrix_ruwen7(ϕ, t, X[1],dP_dT, dP_dTdT,dPhq_dT,dPhq_dalpha, ddPhq_dTdT, ddPhq_dTdalpha,ddPhq_dalphadalpha,n,dn_dmu,dn_dT)

    Ainv = inv(At)

    A_mul_B!(A_i[1], Ainv, Ax)

    jgemvavx!(Source, Ainv, source)



end


function one_d_ideal_matrix_ruwen7(X, tau, r, dP_dT, dP_dTdT,dPhq_dT,dPhq_dalpha, ddPhq_dTdT, ddPhq_dTdalpha,ddPhq_dalphadalpha,n,dn_dmu,dn_dT)

    At = SMatrix{3,3}((X[2] .^2 .*dP_dT  .+ X[2] .^2 .*dPhq_dT  .+ X[1] .*(1  .+ X[2] .^2) .*(dP_dTdT  .+ ddPhq_dTdT),

X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ X[1] .*dP_dTdT  .+ dPhq_dT  .+ X[1] .*ddPhq_dTdT),

sqrt.(1  .+ X[2] .^2) .*dn_dT,

2 .*X[1] .*X[2] .*(dP_dT  .+ dPhq_dT),

(X[1] .*(1  .+ 2 .*X[2] .^2) .*(dP_dT  .+ dPhq_dT)) ./sqrt.(1  .+ X[2] .^2),

(n .*X[2]) ./sqrt.(1  .+ X[2] .^2),

 .-dPhq_dalpha  .+ X[1] .*(1  .+ X[2] .^2) .*ddPhq_dTdalpha,

X[1] .*X[2] .*sqrt.(1  .+ X[2] .^2) .*ddPhq_dTdalpha,

sqrt.(1  .+ X[2] .^2) .*dn_dmu))
    #########################################################################################################################################################################

    Ax = SMatrix{3,3}((X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ X[1] .*dP_dTdT  .+ dPhq_dT  .+ X[1] .*ddPhq_dTdT),

(1  .+ X[2] .^2) .*(dP_dT  .+ dPhq_dT)  .+ X[1] .*X[2] .^2 .*(dP_dTdT  .+ ddPhq_dTdT),

X[2] .*dn_dT,

(X[1] .*(1  .+ 2 .*X[2] .^2) .*(dP_dT  .+ dPhq_dT)) ./sqrt.(1  .+ X[2] .^2),

2 .*X[1] .*X[2] .*(dP_dT  .+ dPhq_dT),

n,

X[1] .*X[2] .*sqrt.(1  .+ X[2] .^2) .*ddPhq_dTdalpha,

dPhq_dalpha  .+ X[1] .*X[2] .^2 .*ddPhq_dTdalpha,

X[2] .*dn_dmu))
    #########################################################################################################################################################################

    source = SVector{3}(((X[1] .*(r  .+ r .*X[2] .^2  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau) .*(dP_dT  .+ dPhq_dT)) ./(r .*tau),

(X[1] .*X[2] .*(r .*sqrt.(1  .+ X[2] .^2)  .+ X[2] .*tau) .*(dP_dT  .+ dPhq_dT)) ./(r .*tau),

n .*(X[2] ./r  .+ sqrt.(1  .+ X[2] .^2) ./tau)))
         return (At,Ax, source)
end

