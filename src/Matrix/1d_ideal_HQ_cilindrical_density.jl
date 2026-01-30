#
# The matrices are based on the notebook: ideal_r_tau_n.nb
# The current is on the EOS of the fluid, i.e. no back reaction of the current, no back feeding of the HQ pressure
# The fields are stored in the order : T, ur, n 
#

function matrix1d_ideal_HQ_density!(A_i, Source, ϕ, t, X, params;free=true)

    #@show t, ϕ
    dP_dT = pressure_derivative(ϕ[1], Val(1), params.eos) #entropy
    dP_dTdT = pressure_derivative(ϕ[1], Val(2), params.eos)


    (At, Ax, source) = one_d_ideal_matrix_ruwen4(ϕ, t, X[1],dP_dT, dP_dTdT)

    Ainv = inv(At)

    A_mul_B!(A_i[1], Ainv, Ax)

    jgemvavx!(Source, Ainv, source)



end


function one_d_ideal_matrix_ruwen4(X, tau, r, dP_dT, dP_dTdT)

    At = SMatrix{3,3}((X[2] .^2 .*dP_dT  .+ X[1] .*(1  .+ X[2] .^2) .*dP_dTdT,

X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ X[1] .*dP_dTdT),

0,

2 .*X[1] .*X[2] .*dP_dT,

(X[1] .*(1  .+ 2 .*X[2] .^2) .*dP_dT) ./sqrt.(1  .+ X[2] .^2),

0,

0,

0,

sqrt.(1  .+ X[2] .^2)))
    #########################################################################################################################################################################

    Ax = SMatrix{3,3}((X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ X[1] .*dP_dTdT),

(1  .+ X[2] .^2) .*dP_dT  .+ X[1] .*X[2] .^2 .*dP_dTdT,

 .-((X[2] .*X[3] .*dP_dTdT) ./( .-(X[2] .^2 .*dP_dT)  .+ X[1] .*(1  .+ X[2] .^2) .*dP_dTdT)),

(X[1] .*(1  .+ 2 .*X[2] .^2) .*dP_dT) ./sqrt.(1  .+ X[2] .^2),

2 .*X[1] .*X[2] .*dP_dT,

(X[1] .*X[3] .*dP_dTdT) ./( .-(X[2] .^2 .*dP_dT)  .+ X[1] .*(1  .+ X[2] .^2) .*dP_dTdT),

0,

0,

X[2]))
    #########################################################################################################################################################################

    source = SVector{3}(((X[1] .*(r  .+ r .*X[2] .^2  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau) .*dP_dT) ./(r .*tau),

(X[1] .*X[2] .*(r .*sqrt.(1  .+ X[2] .^2)  .+ X[2] .*tau) .*dP_dT) ./(r .*tau),

(X[1] .*(1  .+ X[2] .^2) .*X[3] .*(r .*sqrt.(1  .+ X[2] .^2)  .+ X[2] .*tau) .*dP_dTdT) ./(r .*tau .*( .-(X[2] .^2 .*dP_dT)  .+ X[1] .*(1  .+ X[2] .^2) .*dP_dTdT))))
         return (At,Ax, source)
end

