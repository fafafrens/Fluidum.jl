#
# The matrices are based on the notebook: ideal_r_tau_br_p_n.nb
# The fields are stored in the order : T, ur, n 
# HQ back reaction and pressure feedback included
#

function matrix1d_ideal_HQ_br_p_density!(A_i, Source, ϕ, t, X, params;free=true)

    #@show t, ϕ
    dP_dT = pressure_derivative(ϕ[1], Val(1), params.eos) #entropy
    dP_dTdT = pressure_derivative(ϕ[1], Val(2), params.eos)



    if free == true 
        α = invert_n_for_alpha(ϕ[1], ϕ[3]; m=params.diffusion.mass)
        therm = hq_pressure(ϕ[1], α; m=params.diffusion.mass)    
        dPhq_dT, dPhq_dalpha = therm.pressure_derivative
        ddPhq_dTdT, ddPhq_dTdalpha, ddPhq_dalphadalpha = therm.pressure_hessian

        #therm ,_ = hq_pressure_Tn(ϕ[1], ϕ[3]; m=params.diffusion.mass)    
        #dPhq_dT, dPhq_dn = therm.pressure_derivative
        #ddPhq_dTdT, ddPhq_dTdn, ddPhq_dndn = therm.pressure_hessian
    else 
       @warn "not free ideal not implemented."
    end

    (At, Ax, source) = one_d_ideal_matrix_ruwen6(ϕ, t, X[1],dP_dT, dP_dTdT,dPhq_dT,dPhq_dalpha, ddPhq_dTdT, ddPhq_dTdalpha,ddPhq_dalphadalpha)

    Ainv = inv(At)

    A_mul_B!(A_i[1], Ainv, Ax)

    jgemvavx!(Source, Ainv, source)



end


function one_d_ideal_matrix_ruwen6(X, tau, r, dP_dT, dP_dTdT,dPhq_dT,dPhq_dalpha, ddPhq_dTdT, ddPhq_dTdalpha,ddPhq_dalphadalpha)

    At = SMatrix{3,3}((( .-dPhq_dalpha .^2  .+ X[1] .*(2  .+ X[2] .^2) .*dPhq_dalpha .*ddPhq_dTdalpha  .+ X[1] .*(X[2] .^2 .*dP_dT .*ddPhq_dalphadalpha  .+ X[1] .*(1  .+ X[2] .^2) .*dP_dTdT .*ddPhq_dalphadalpha  .+ X[2] .^2 .*ddPhq_dalphadalpha .*dPhq_dT  .- X[1] .*ddPhq_dTdalpha .^2  .- X[1] .*X[2] .^2 .*ddPhq_dTdalpha .^2  .+ X[1] .*ddPhq_dalphadalpha .*ddPhq_dTdT  .+ X[1] .*X[2] .^2 .*ddPhq_dalphadalpha .*ddPhq_dTdT)) ./(X[1] .*ddPhq_dalphadalpha),

(X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT .*ddPhq_dalphadalpha  .+ X[1] .*dP_dTdT .*ddPhq_dalphadalpha  .+ ddPhq_dalphadalpha .*dPhq_dT  .+ dPhq_dalpha .*ddPhq_dTdalpha  .- X[1] .*ddPhq_dTdalpha .^2  .+ X[1] .*ddPhq_dalphadalpha .*ddPhq_dTdT)) ./ddPhq_dalphadalpha,

0,

2 .*X[1] .*X[2] .*(dP_dT  .+ dPhq_dT),

(X[1] .*(1  .+ 2 .*X[2] .^2) .*(dP_dT  .+ dPhq_dT)) ./sqrt.(1  .+ X[2] .^2),

(X[2] .*X[3]) ./sqrt.(1  .+ X[2] .^2),

(X[1] .^2 .*X[2] .*sqrt.(1  .+ X[2] .^2) .*ddPhq_dTdalpha) ./ddPhq_dalphadalpha,

(X[1] .*(dPhq_dalpha  .+ X[1] .*X[2] .^2 .*ddPhq_dTdalpha)) ./ddPhq_dalphadalpha,

sqrt.(1  .+ X[2] .^2)))
    #########################################################################################################################################################################

    Ax = SMatrix{3,3}(((X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT .*ddPhq_dalphadalpha  .+ X[1] .*dP_dTdT .*ddPhq_dalphadalpha  .+ ddPhq_dalphadalpha .*dPhq_dT  .+ dPhq_dalpha .*ddPhq_dTdalpha  .- X[1] .*ddPhq_dTdalpha .^2  .+ X[1] .*ddPhq_dalphadalpha .*ddPhq_dTdT)) ./ddPhq_dalphadalpha,

(dPhq_dalpha .^2  .+ X[1] .*( .-1  .+ X[2] .^2) .*dPhq_dalpha .*ddPhq_dTdalpha  .+ X[1] .*((1  .+ X[2] .^2) .*dP_dT .*ddPhq_dalphadalpha  .+ X[1] .*X[2] .^2 .*dP_dTdT .*ddPhq_dalphadalpha  .+ ddPhq_dalphadalpha .*dPhq_dT  .+ X[2] .^2 .*ddPhq_dalphadalpha .*dPhq_dT  .- X[1] .*X[2] .^2 .*ddPhq_dTdalpha .^2  .+ X[1] .*X[2] .^2 .*ddPhq_dalphadalpha .*ddPhq_dTdT)) ./(X[1] .*ddPhq_dalphadalpha),

0,

(X[1] .*(1  .+ 2 .*X[2] .^2) .*(dP_dT  .+ dPhq_dT)) ./sqrt.(1  .+ X[2] .^2),

2 .*X[1] .*X[2] .*(dP_dT  .+ dPhq_dT),

X[3],

(X[1] .*( .-dPhq_dalpha  .+ X[1] .*(1  .+ X[2] .^2) .*ddPhq_dTdalpha)) ./ddPhq_dalphadalpha,

(X[1] .^2 .*X[2] .*sqrt.(1  .+ X[2] .^2) .*ddPhq_dTdalpha) ./ddPhq_dalphadalpha,

X[2]))
    #########################################################################################################################################################################

    source = SVector{3}(((X[1] .*(r  .+ r .*X[2] .^2  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau) .*(dP_dT  .+ dPhq_dT)) ./(r .*tau),

(X[1] .*X[2] .*(r .*sqrt.(1  .+ X[2] .^2)  .+ X[2] .*tau) .*(dP_dT  .+ dPhq_dT)) ./(r .*tau),

X[3] .*(X[2] ./r  .+ sqrt.(1  .+ X[2] .^2) ./tau)))
         return (At,Ax, source)
end

