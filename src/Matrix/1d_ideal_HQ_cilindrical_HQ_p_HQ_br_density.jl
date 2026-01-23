#
# The matrices are based on the notebook: O(2)MIS_viscious_r_tau_all_ideal_br_p_n.nb
# The fields are stored in the order : T, ur, n 
#
function matrix1d_full_ideal_HQ_br_HQ_p_density!(A_i, Source, ϕ, t, X, params;free=true)

    #@show t, ϕ
    dP_dT = pressure_derivative(ϕ[1], Val(1), params.eos) #entropy
    dP_dTdT = pressure_derivative(ϕ[1], Val(2), params.eos)



    if free == true 
        therm ,_ = hq_pressure_Tn(ϕ[1], ϕ[3]; m=params.diffusion.mass)    
        dPhq_dT, dPhq_dn = therm.pressure_derivative
        ddPhq_dTdT, ddPhq_dTdn, ddPhq_dndn = therm.pressure_hessian
    else 
       @warn "not free ideal not implemented."
    end

    (At, Ax, source) = one_d_ideal_matrix_ruwen3(ϕ, t, X[1],dP_dT, dP_dTdT,dPhq_dT,dPhq_dn, ddPhq_dTdT, ddPhq_dTdn)

    Ainv = inv(At)

    A_mul_B!(A_i[1], Ainv, Ax)

    jgemvavx!(Source, Ainv, source)



end


function one_d_ideal_matrix_ruwen3(X, tau, r, dP_dT, dP_dTdT,dPhq_dT,dPhq_dn, ddPhq_dTdT, ddPhq_dTdn)

    At = SMatrix{3,3}((X[2] .^2 .*dP_dT  .+ X[2] .^2 .*dPhq_dT  .+ X[1] .*(1  .+ X[2] .^2) .*(dP_dTdT  .+ ddPhq_dTdT),

X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ X[1] .*dP_dTdT  .+ dPhq_dT  .+ X[1] .*ddPhq_dTdT),

0,

2 .*X[1] .*X[2] .*(dP_dT  .+ dPhq_dT),

(X[1] .*(1  .+ 2 .*X[2] .^2) .*(dP_dT  .+ dPhq_dT)) ./sqrt.(1  .+ X[2] .^2),

(X[2] .*X[3]) ./sqrt.(1  .+ X[2] .^2),

 .-dPhq_dn  .+ X[1] .*(1  .+ X[2] .^2) .*ddPhq_dTdn,

X[1] .*X[2] .*sqrt.(1  .+ X[2] .^2) .*ddPhq_dTdn,

sqrt.(1  .+ X[2] .^2)))
    #########################################################################################################################################################################

    Ax = SMatrix{3,3}((X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ X[1] .*dP_dTdT  .+ dPhq_dT  .+ X[1] .*ddPhq_dTdT),

(1  .+ X[2] .^2) .*(dP_dT  .+ dPhq_dT)  .+ X[1] .*X[2] .^2 .*(dP_dTdT  .+ ddPhq_dTdT),

0,

(X[1] .*(1  .+ 2 .*X[2] .^2) .*(dP_dT  .+ dPhq_dT)) ./sqrt.(1  .+ X[2] .^2),

2 .*X[1] .*X[2] .*(dP_dT  .+ dPhq_dT),

X[3],

X[1] .*X[2] .*sqrt.(1  .+ X[2] .^2) .*ddPhq_dTdn,

dPhq_dn  .+ X[1] .*X[2] .^2 .*ddPhq_dTdn,

X[2]))
    #########################################################################################################################################################################

    source = SVector{3}(((X[1] .*(r  .+ r .*X[2] .^2  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau) .*(dP_dT  .+ dPhq_dT)) ./(r .*tau),

(X[1] .*X[2] .*(r .*sqrt.(1  .+ X[2] .^2)  .+ X[2] .*tau) .*(dP_dT  .+ dPhq_dT)) ./(r .*tau),

X[3] .*(X[2] ./r  .+ sqrt.(1  .+ X[2] .^2) ./tau)))
         return (At,Ax, source)
end

