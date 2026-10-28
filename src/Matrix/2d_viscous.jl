# the convention here are T, ux, uy, \[Pi]yy, \[Pi]zz, \[Pi]xy, \[Pi]B
# two_d_viscous_matrix_derived (end of this file, derived symbolically) is the default; the legacy
# two_d_viscous_matrix below has ten wrong entries in the velocity columns of the pi^yy and pi^xy rows (legacy = true)

"""
    matrix2d_visc!(; legacy = false, delta_pipi = true)
2+1D viscous hydro, derived matrix; `delta_pipi = false` drops the DNMR term, `legacy = true` uses the legacy matrix.
"""
function matrix2d_visc!(; legacy::Bool = false, delta_pipi::Bool = true)
    (A_i, Source, ϕ, t, X, params) -> matrix2d_visc!(A_i, Source, ϕ, t, X, params;
                                                     legacy = legacy, delta_pipi = delta_pipi)
end

@inbounds @fastmath function matrix2d_visc!(A_i,Source,ϕ,t,X,params; legacy::Bool = false, delta_pipi::Bool = true)
 

    therm=thermodynamic(ϕ[1],params)
        
    etaVisc=viscosity(ϕ[1],therm,params.shear)
    
    tauS=τ_shear(ϕ[1],therm,params.shear)
    
    tauB=τ_bulk(ϕ[1],therm,params.bulk)
    
    zeta=bulk_viscosity(ϕ[1],therm,params.bulk)
    
 
    (At,Ax, Ay, source)= legacy ? two_d_viscous_matrix(ϕ,t,therm.pressure,therm.pressure_derivative[1],therm.pressure_hessian[1],zeta,etaVisc,tauS,tauB) :
                       delta_pipi ? two_d_viscous_matrix_derived(ϕ,t,therm.pressure,therm.pressure_derivative[1],therm.pressure_hessian[1],zeta,etaVisc,tauS,tauB) :
                                    two_d_viscous_matrix_derived(ϕ,t,therm.pressure,therm.pressure_derivative[1],therm.pressure_hessian[1],zeta,etaVisc,tauS,tauB,zero(tauS))
    

        
        Ainv= inv(At)

    
        A_mul_B!(A_i[1], Ainv,Ax)
        A_mul_B!(A_i[2], Ainv,Ay)
        
        jgemvavx!(Source, Ainv,source)
    
        
    
    end 



@inbounds @fastmath function two_d_viscous_matrix(u,tau,p,dtp,dtdtp,zeta,visc,tauS,tauB)
    
At=SMatrix{7,7}(

    dtp*(^(u[2],2) + ^(u[3],2)) + dtdtp*u[1]*(1 + ^(u[2],2) + ^(u[3],2))

    ,(dtp + dtdtp*u[1])*u[2]*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,(dtp + dtdtp*u[1])*u[3]*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,0

    ,0

    ,0

    ,0

    ,(2*(dtp*u[1]*u[2]*(1 + ^(u[3],2)) + u[3]*u[6] + u[2]*(u[7] + ^(u[3],2)*u[7] - u[4] - u[5])))/(1 + ^(u[3],2))

    ,(dtp*u[1]*(1 + ^(u[3],2))*(1 + ^(u[2],2) + ^(u[3],2))*(1 + 2*^(u[2],2) + ^(u[3],2)) + (1 + ^(u[3],2))*(1 + ^(u[2],2) + ^(u[3],2))*(1 + 2*^(u[2],2) + ^(u[3],2))*u[7] + 3*u[2]*u[3]*u[6] + 2*^(u[2],3)*u[3]*u[6] + 3*u[2]*^(u[3],3)*u[6] - u[4] - 3*^(u[2],2)*u[4] - 2*^(u[2],4)*u[4] - ^(u[3],2)*u[4] - 3*^(u[2],2)*^(u[3],2)*u[4] - (1 + ^(u[2],2) + ^(u[3],2))*(1 + 2*^(u[2],2) + ^(u[3],2))*u[5])/((1 + ^(u[3],2))*^(1 + ^(u[2],2) + ^(u[3],2),1.5))

    ,(u[6] + u[3]*(dtp*u[1]*u[2]*(1 + ^(u[2],2) + ^(u[3],2)) + ^(u[2],3)*u[7] + u[3]*u[6] + u[2]*(u[7] + ^(u[3],2)*u[7] - u[4])))/^(1 + ^(u[2],2) + ^(u[3],2),1.5)

    ,(-2*(1 + ^(u[3],2))*(u[2]*visc + 3*tauS*u[3]*u[6]) + 6*tauS*u[2]*^(u[3],2)*u[4])/(3*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,(-2*u[2]*visc)/(3*^(tau,2)*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,(-3*u[2]*(-1 + tauS*(2 + 4*^(u[3],2)))*u[6] + 3*(1 + 2*tauS)*u[3]*u[4] + 6*tauS*(u[3] + ^(u[3],3))*u[5] + 2*^(u[2],2)*u[3]*(-2*visc + 3*tauS*(2*u[4] + u[5])))/(6*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,(u[2]*zeta)/sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,(2*(u[2]*u[6] + u[3]*(dtp*u[1]*^(1 + ^(u[3],2),2) + ^(1 + ^(u[3],2),2)*u[7] + u[4] + u[2]*(-(u[3]*u[6]) + u[2]*(u[4] + u[5])))))/^(1 + ^(u[3],2),2)

    ,(dtp*u[1]*u[2]*u[3]*^(1 + ^(u[3],2),2)*(1 + ^(u[2],2) + ^(u[3],2)) - 2*^(u[2],4)*(-1 + ^(u[3],2))*u[6] + ^(1 + ^(u[3],2),2)*u[6] - 3*^(u[2],2)*(-1 + ^(u[3],4))*u[6] + 2*^(u[2],5)*u[3]*(u[4] + u[5]) + u[2]*u[3]*(1 + ^(u[3],2))*(^(1 + ^(u[3],2),2)*u[7] + 3*u[4] + u[5] + ^(u[3],2)*u[5]) + ^(u[2],3)*u[3]*(^(1 + ^(u[3],2),2)*u[7] + 5*u[4] + 3*u[5] + 3*^(u[3],2)*(u[4] + u[5])))/(^(1 + ^(u[3],2),2)*^(1 + ^(u[2],2) + ^(u[3],2),1.5))

    ,(dtp*u[1]*(1 + ^(u[2],2) + ^(u[3],2))*(1 + ^(u[2],2) + 2*^(u[3],2)) + (1 + ^(u[2],2) + ^(u[3],2))*(1 + ^(u[2],2) + 2*^(u[3],2))*u[7] - u[2]*u[3]*u[6] + u[4] + ^(u[2],2)*u[4])/^(1 + ^(u[2],2) + ^(u[3],2),1.5)

    ,(4*^(u[3],3)*visc + 3*u[2]*u[6] + 6*tauS*u[2]*^(u[3],2)*u[6] + u[3]*((4 + 6*^(u[2],2))*visc + 3*u[4] - 6*tauS*(1 + ^(u[2],2))*u[4]))/(3*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,(-2*u[3]*visc)/(3*^(tau,2)*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,(-3*(-1 + 2*tauS)*u[3]*(1 + ^(u[3],2))*u[6] + 6*^(u[2],2)*u[3]*(1 + 2*tauS*^(u[3],2))*u[6] - 3*^(u[2],3)*(-4*(1 + ^(u[3],2))*visc + u[4] + 2*tauS*u[4] + u[5] + 2*tauS*^(u[3],2)*(2*u[4] + u[5])) + u[2]*(4*(3 + 5*^(u[3],2) + 2*^(u[3],4))*visc - 3*(u[4] + 2*tauS*u[4] + u[5]) - 3*^(u[3],2)*(u[5] + 2*tauS*(2*u[4] + u[5] + ^(u[3],2)*u[5]))))/(6*(1 + ^(u[3],2))*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,(u[3]*zeta)/sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,(-^(u[2],2) + ^(u[3],2))/(1 + ^(u[3],2))

    ,-((u[2]*(1 + ^(u[2],2)))/((1 + ^(u[3],2))*sqrt(1 + ^(u[2],2) + ^(u[3],2))))

    ,u[3]/sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,tauS*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,0

    ,0

    ,0

    ,-(^(u[2],2)/(1 + ^(u[3],2)))

    ,-((u[2]*sqrt(1 + ^(u[2],2) + ^(u[3],2)))/(1 + ^(u[3],2)))

    ,0

    ,0

    ,(tauS*sqrt(1 + ^(u[2],2) + ^(u[3],2)))/^(tau,2)

    ,0

    ,0

    ,(2*u[2]*u[3])/(1 + ^(u[3],2))

    ,(u[3] + 2*^(u[2],2)*u[3] + ^(u[3],3))/((1 + ^(u[3],2))*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,u[2]/sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,0

    ,0

    ,tauS*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,0

    ,^(u[2],2) + ^(u[3],2)

    ,u[2]*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,u[3]*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,0

    ,0

    ,0

    ,tauB*sqrt(1 + ^(u[2],2) + ^(u[3],2))

)
    #########################################################################################################################################################################
    
Ax=SMatrix{7,7}(

    (dtp + dtdtp*u[1])*u[2]*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,dtp + (dtp + dtdtp*u[1])*^(u[2],2)

    ,(dtp + dtdtp*u[1])*u[2]*u[3]

    ,0

    ,0

    ,0

    ,0

    ,(dtp*u[1]*(1 + ^(u[3],2))*(1 + ^(u[2],2) + ^(u[3],2))*(1 + 2*^(u[2],2) + ^(u[3],2)) + (1 + ^(u[3],2))*(1 + ^(u[2],2) + ^(u[3],2))*(1 + 2*^(u[2],2) + ^(u[3],2))*u[7] + 3*u[2]*u[3]*u[6] + 2*^(u[2],3)*u[3]*u[6] + 3*u[2]*^(u[3],3)*u[6] - u[4] - 3*^(u[2],2)*u[4] - 2*^(u[2],4)*u[4] - ^(u[3],2)*u[4] - 3*^(u[2],2)*^(u[3],2)*u[4] - (1 + ^(u[2],2) + ^(u[3],2))*(1 + 2*^(u[2],2) + ^(u[3],2))*u[5])/((1 + ^(u[3],2))*^(1 + ^(u[2],2) + ^(u[3],2),1.5))

    ,(2*(dtp*u[1]*u[2]*(1 + ^(u[3],2)) + u[3]*u[6] + u[2]*(u[7] + ^(u[3],2)*u[7] - u[4] - u[5])))/(1 + ^(u[3],2))

    ,u[3]*(dtp*u[1] + u[7])

    ,(-2*(1 + ^(u[3],2))*visc)/3. + (2*tauS*u[2]*u[3]*(-((1 + ^(u[3],2))*u[6]) + u[2]*u[3]*u[4]))/(1 + ^(u[2],2) + ^(u[3],2))

    ,(-2*visc)/(3*^(tau,2))

    ,(u[2]*(3*(1 - 2*tauS)*u[2]*u[6] - 12*tauS*u[2]*^(u[3],2)*u[6] + ^(u[3],3)*(-4*visc + 6*tauS*u[5]) + u[3]*(-4*(1 + ^(u[2],2))*visc + 3*(1 + 2*tauS)*u[4] + 6*tauS*(u[5] + ^(u[2],2)*(2*u[4] + u[5])))))/(6*(1 + ^(u[2],2) + ^(u[3],2)))

    ,zeta

    ,(dtp*u[1]*u[2]*u[3]*^(1 + ^(u[3],2),2)*(1 + ^(u[2],2) + ^(u[3],2)) - 2*^(u[2],4)*(-1 + ^(u[3],2))*u[6] + ^(1 + ^(u[3],2),2)*u[6] - 3*^(u[2],2)*(-1 + ^(u[3],4))*u[6] + 2*^(u[2],5)*u[3]*(u[4] + u[5]) + u[2]*u[3]*(1 + ^(u[3],2))*(^(1 + ^(u[3],2),2)*u[7] + 3*u[4] + u[5] + ^(u[3],2)*u[5]) + ^(u[2],3)*u[3]*(^(1 + ^(u[3],2),2)*u[7] + 5*u[4] + 3*u[5] + 3*^(u[3],2)*(u[4] + u[5])))/(^(1 + ^(u[3],2),2)*^(1 + ^(u[2],2) + ^(u[3],2),1.5))

    ,(2*(u[2]*(u[6] - ^(u[3],2)*u[6]) + u[3]*u[4] + ^(u[2],2)*u[3]*(u[4] + u[5])))/^(1 + ^(u[3],2),2)

    ,u[2]*(dtp*u[1] + u[7])

    ,u[6] + (2*u[2]*u[3]*((1 + ^(u[2],2) + ^(u[3],2))*visc + tauS*u[2]*u[3]*u[6]) - 2*tauS*(u[2] + ^(u[2],3))*u[3]*u[4])/(1 + ^(u[2],2) + ^(u[3],2))

    ,0

    ,-(-4*(1 + ^(u[2],2))*(1 + ^(u[3],2))*(1 + ^(u[2],2) + ^(u[3],2))*visc + (-3 + 2*tauS)*u[2]*u[3]*(1 + ^(u[3],2))*u[6] - 2*^(u[2],3)*u[3]*(1 + 2*tauS*^(u[3],2))*u[6] + (1 + ^(u[3],2))*(2*u[4] + u[5] + ^(u[3],2)*u[5]) + ^(u[2],4)*(u[4] + tauS*(2 + 4*^(u[3],2))*u[4] + u[5] + 2*tauS*^(u[3],2)*u[5]) + ^(u[2],2)*((3 + 2*tauS + 2*(1 + 2*tauS)*^(u[3],2))*u[4] + 2*(1 + ^(u[3],2))*(1 + tauS*^(u[3],2))*u[5]))/(2*(1 + ^(u[3],2))*(1 + ^(u[2],2) + ^(u[3],2)))

    ,0

    ,-((u[2]*(1 + ^(u[2],2)))/((1 + ^(u[3],2))*sqrt(1 + ^(u[2],2) + ^(u[3],2))))

    ,-((1 + ^(u[2],2))/(1 + ^(u[3],2)))

    ,0

    ,tauS*u[2]

    ,0

    ,0

    ,0

    ,-((u[2]*sqrt(1 + ^(u[2],2) + ^(u[3],2)))/(1 + ^(u[3],2)))

    ,-1 - ^(u[2],2)/(1 + ^(u[3],2))

    ,0

    ,0

    ,(tauS*u[2])/^(tau,2)

    ,0

    ,0

    ,(u[3] + 2*^(u[2],2)*u[3] + ^(u[3],3))/((1 + ^(u[3],2))*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,(2*u[2]*u[3])/(1 + ^(u[3],2))

    ,1

    ,0

    ,0

    ,tauS*u[2]

    ,0

    ,u[2]*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,1 + ^(u[2],2)

    ,u[2]*u[3]

    ,0

    ,0

    ,0

    ,tauB*u[2]

)
    #########################################################################################################################################################################
Ay=SMatrix{7,7}(

    (dtp + dtdtp*u[1])*u[3]*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,(dtp + dtdtp*u[1])*u[2]*u[3]

    ,dtp + (dtp + dtdtp*u[1])*^(u[3],2)

    ,0

    ,0

    ,0

    ,0

    ,(u[6] + u[3]*(dtp*u[1]*u[2]*(1 + ^(u[2],2) + ^(u[3],2)) + ^(u[2],3)*u[7] + u[3]*u[6] + u[2]*(u[7] + ^(u[3],2)*u[7] - u[4])))/^(1 + ^(u[2],2) + ^(u[3],2),1.5)

    ,u[3]*(dtp*u[1] + u[7])

    ,0

    ,-(((1 + 2*tauS*^(u[3],2))*(u[6] + ^(u[3],2)*u[6] - u[2]*u[3]*u[4]))/(1 + ^(u[2],2) + ^(u[3],2)))

    ,0

    ,(-(u[2]*u[3]*(1 + tauS*(2 + 4*^(u[3],2)))*u[6]) + 2*u[4] + u[5] + ^(u[2],2)*(1 + 2*tauS*^(u[3],2))*(2*u[4] + u[5]) + ^(u[3],2)*(u[4] + 2*tauS*u[4] + u[5] + 2*tauS*(1 + ^(u[3],2))*u[5]))/(2*(1 + ^(u[2],2) + ^(u[3],2)))

    ,0

    ,(dtp*u[1]*(1 + ^(u[2],2) + ^(u[3],2))*(1 + ^(u[2],2) + 2*^(u[3],2)) + (1 + ^(u[2],2) + ^(u[3],2))*(1 + ^(u[2],2) + 2*^(u[3],2))*u[7] - u[2]*u[3]*u[6] + u[4] + ^(u[2],2)*u[4])/^(1 + ^(u[2],2) + ^(u[3],2),1.5)

    ,u[2]*(dtp*u[1] + u[7])

    ,2*u[3]*(dtp*u[1] + u[7])

    ,(4*(1 + ^(u[3],2))*visc)/3. + (u[3]*(u[2]*(1 + 2*tauS*^(u[3],2))*u[6] + (1 - 2*tauS*(1 + ^(u[2],2)))*u[3]*u[4]))/(1 + ^(u[2],2) + ^(u[3],2))

    ,(-2*visc)/(3*^(tau,2))

    ,(u[3]*(-3*(-1 + 2*tauS)*u[3]*(1 + ^(u[3],2))*u[6] + 6*^(u[2],2)*u[3]*(1 + 2*tauS*^(u[3],2))*u[6] + ^(u[2],3)*(8*(1 + ^(u[3],2))*visc - 6*tauS*^(u[3],2)*(2*u[4] + u[5]) - 3*(u[4] + 2*tauS*u[4] + u[5])) + u[2]*(8*^(1 + ^(u[3],2),2)*visc - 3*(u[4] + tauS*(2 + 4*^(u[3],2))*u[4] + (1 + ^(u[3],2))*(1 + 2*tauS*^(u[3],2))*u[5]))))/(6*(1 + ^(u[3],2))*(1 + ^(u[2],2) + ^(u[3],2)))

    ,zeta

    ,u[3]/sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,0

    ,1

    ,tauS*u[3]

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,(tauS*u[3])/^(tau,2)

    ,0

    ,0

    ,u[2]/sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,1

    ,0

    ,0

    ,0

    ,tauS*u[3]

    ,0

    ,u[3]*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,u[2]*u[3]

    ,1 + ^(u[3],2)

    ,0

    ,0

    ,0

    ,tauB*u[3]

)
    #########################################################################################################################################################


    #########################################################################################################################################################

source=SVector{7}(

    (dtp*u[1]*(1 + ^(u[3],2))*(1 + ^(u[2],2) + ^(u[3],2)) + u[7] + 2*u[2]*u[3]*u[6] + ^(u[2],2)*(u[7] + ^(u[3],2)*u[7] - u[4] - u[5]) + u[5] + ^(u[3],2)*((2 + ^(u[3],2))*u[7] + u[4] + u[5]))/(tau*(1 + ^(u[3],2)))

    ,(dtp*u[1]*u[2]*(1 + ^(u[3],2))*(1 + ^(u[2],2) + ^(u[3],2)) + 2*^(u[2],2)*u[3]*u[6] + u[3]*(1 + ^(u[3],2))*u[6] + ^(u[2],3)*(u[7] + ^(u[3],2)*u[7] - u[4] - u[5]) + u[2]*(^(1 + ^(u[3],2),2)*u[7] - u[4] - (1 + ^(u[3],2))*u[5]))/(tau*(1 + ^(u[3],2))*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,(dtp*u[1]*u[3]*(1 + ^(u[2],2) + ^(u[3],2)) + ^(u[3],3)*u[7] + u[2]*u[6] + u[3]*(u[7] + ^(u[2],2)*u[7] + u[4]))/(tau*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,(-2*(1 + ^(u[3],2))*sqrt(1 + ^(u[2],2) + ^(u[3],2))*visc)/(3*tau) + u[4]

    ,(4*sqrt(1 + ^(u[2],2) + ^(u[3],2))*visc + 3*tau*u[5])/(3*^(tau,3))

    ,(-2*u[2]*u[3]*sqrt(1 + ^(u[2],2) + ^(u[3],2))*visc)/(3*tau) + u[6]

    ,(sqrt(1 + ^(u[2],2) + ^(u[3],2))*zeta)/tau + u[7]

)

return (At,Ax,Ay,source)
end


# generated by symbolic computation (Mathematica)
"""
    two_d_viscous_matrix_derived(u, tau, p, dtp, dtdtp, zeta, visc, tauS, tauB, dpp = delta_pipi(tauS))
(A_t, A_x, A_y, src) for phi = (T, u^x, u^y, pi^{yy}, pi^eta_eta, pi^{xy}, Pi), the layout of `two_d_viscous_matrix`.
Israel-Stewart + the DNMR term dpp pi theta, source on the LHS; dpp = 0 is plain Israel-Stewart.
"""
function two_d_viscous_matrix_derived(u, tau, pres, dtp, dtdtp, zeta, visc, tauS, tauB, dpp = delta_pipi(tauS))
    T, ux, uy, piyy, pizz, pixy, piB = u[1], u[2], u[3], u[4], u[5], u[6], u[7]
    At = SMatrix{7,7}((
        -dtp + (dtp + dtdtp*T)*(1 + ux^2 + uy^2),
        (dtp + dtdtp*T)*ux*sqrt(1 + ux^2 + uy^2),
        (dtp + dtdtp*T)*uy*sqrt(1 + ux^2 + uy^2),
        0,
        0,
        0,
        0,
        (2*(-(pizz*ux) - piyy*ux + dtp*T*ux + pixy*uy + dtp*T*ux*uy^2 + piB*ux*(1 + uy^2)))/(1 + uy^2),
        ((piB + dtp*T)*ux^2*(1 + ux^2 + uy^2) + (piB + dtp*T)*(1 + ux^2 + uy^2)^2 + (ux*(piyy*(ux + ux^3) + pizz*ux*(1 + ux^2 + uy^2) - pixy*uy*(1 + 2*ux^2 + uy^2)))/(1 + uy^2) - ((1 + ux^2 + uy^2)*(piyy + 3*piyy*ux^2 - 4*pixy*ux*uy + pizz*(1 + 3*ux^2 + uy^2)))/(1 + uy^2))/(1 + ux^2 + uy^2)^(3/2),
        (pixy*(1 + uy^2) + ux*uy*(-piyy + piB*(1 + ux^2 + uy^2) + dtp*T*(1 + ux^2 + uy^2)))/(1 + ux^2 + uy^2)^(3/2),
        (-2*(3*pixy*tauS*uy*(1 + uy^2) + ux*(visc + uy^2*(-3*piyy*tauS + visc))))/(3*sqrt(1 + ux^2 + uy^2)),
        (-2*ux*visc)/(3*tau^2*sqrt(1 + ux^2 + uy^2)),
        (-3*pixy*tauS*ux*(1 + 2*uy^2) + uy*(3*piyy*(tauS + 2*tauS*ux^2) + 3*pizz*tauS*(1 + ux^2 + uy^2) + (3 + ux^2 + 3*uy^2)*visc))/(3*sqrt(1 + ux^2 + uy^2)),
        (ux*zeta)/sqrt(1 + ux^2 + uy^2),
        (2*(pixy*(ux - ux*uy^2) + uy*(piyy + dtp*T + pizz*ux^2 + piyy*ux^2 + 2*dtp*T*uy^2 + dtp*T*uy^4 + piB*(1 + uy^2)^2)))/(1 + uy^2)^2,
        ((piB + dtp*T)*ux*uy*(1 + ux^2 + uy^2) + ((1 + ux^2 + uy^2)*(2*ux*(piyy + pizz*ux^2 + piyy*ux^2)*uy + pixy*(-2*ux^2*(-1 + uy^2) + (1 + uy^2)^2)))/(1 + uy^2)^2 + uy*(-(pixy*uy) + (ux*(piyy + piyy*ux^2 - 2*pixy*ux*uy + pizz*(1 + ux^2 + uy^2)))/(1 + uy^2)))/(1 + ux^2 + uy^2)^(3/2),
        (-(uy*(pixy*ux + piyy*uy)) + piyy*(1 + ux^2 + uy^2) + (piB + dtp*T)*uy^2*(1 + ux^2 + uy^2) + (piB + dtp*T)*(1 + ux^2 + uy^2)^2)/(1 + ux^2 + uy^2)^(3/2),
        (2*uy*(-3*piyy*tauS*(1 + ux^2) + 3*pixy*tauS*ux*uy + (2 + 3*ux^2 + 2*uy^2)*visc))/(3*sqrt(1 + ux^2 + uy^2)),
        (-2*uy*visc)/(3*tau^2*sqrt(1 + ux^2 + uy^2)),
        (-3*piyy*tauS*ux*(1 + ux^2)*(1 + 2*uy^2) + 3*pixy*tauS*uy*(-1 + (-1 + 2*ux^2)*uy^2) + ux*(-3*pizz*tauS*uy^2*(1 + ux^2 + uy^2) + (1 + uy^2)*(3 + 3*ux^2 + uy^2)*visc))/(3*(1 + uy^2)*sqrt(1 + ux^2 + uy^2)),
        (uy*zeta)/sqrt(1 + ux^2 + uy^2),
        (-ux^2 + uy^2)/(1 + uy^2),
        -((ux*(1 + ux^2))/((1 + uy^2)*sqrt(1 + ux^2 + uy^2))),
        uy/sqrt(1 + ux^2 + uy^2),
        tauS*sqrt(1 + ux^2 + uy^2),
        0,
        0,
        0,
        -(ux^2/(1 + uy^2)),
        -((ux*sqrt(1 + ux^2 + uy^2))/(1 + uy^2)),
        0,
        0,
        (tauS*sqrt(1 + ux^2 + uy^2))/tau^2,
        0,
        0,
        (2*ux*uy)/(1 + uy^2),
        (uy*(1 + (2*ux^2)/(1 + uy^2)))/sqrt(1 + ux^2 + uy^2),
        ux/sqrt(1 + ux^2 + uy^2),
        0,
        0,
        tauS*sqrt(1 + ux^2 + uy^2),
        0,
        ux^2 + uy^2,
        ux*sqrt(1 + ux^2 + uy^2),
        uy*sqrt(1 + ux^2 + uy^2),
        0,
        0,
        0,
        tauB*sqrt(1 + ux^2 + uy^2)))

    Ax = SMatrix{7,7}((
        (dtp + dtdtp*T)*ux*sqrt(1 + ux^2 + uy^2),
        dtp + dtp*ux^2 + dtdtp*T*ux^2,
        (dtp + dtdtp*T)*ux*uy,
        0,
        0,
        0,
        0,
        ((piB + dtp*T)*ux^2*(1 + ux^2 + uy^2) + (piB + dtp*T)*(1 + ux^2 + uy^2)^2 + (ux*(piyy*(ux + ux^3) + pizz*ux*(1 + ux^2 + uy^2) - pixy*uy*(1 + 2*ux^2 + uy^2)))/(1 + uy^2) - ((1 + ux^2 + uy^2)*(piyy + 3*piyy*ux^2 - 4*pixy*ux*uy + pizz*(1 + 3*ux^2 + uy^2)))/(1 + uy^2))/(1 + ux^2 + uy^2)^(3/2),
        2*(piB + dtp*T)*ux - (2*(pizz*ux + piyy*ux - pixy*uy))/(1 + uy^2),
        (piB + dtp*T)*uy,
        (-2*(-3*piyy*tauS*ux^2*uy^2 + 3*pixy*tauS*ux*uy*(1 + uy^2) + (1 + uy^2)*(1 + ux^2 + uy^2)*visc))/(3*(1 + ux^2 + uy^2)),
        (-2*visc)/(3*tau^2),
        (ux*(-3*pixy*tauS*ux*(1 + 2*uy^2) + uy*(3*piyy*(tauS + 2*tauS*ux^2) + 3*pizz*tauS*(1 + ux^2 + uy^2) + (1 + ux^2 + uy^2)*visc)))/(3*(1 + ux^2 + uy^2)),
        zeta,
        ((piB + dtp*T)*ux*uy*(1 + ux^2 + uy^2) + ((1 + ux^2 + uy^2)*(2*ux*(piyy + pizz*ux^2 + piyy*ux^2)*uy + pixy*(-2*ux^2*(-1 + uy^2) + (1 + uy^2)^2)))/(1 + uy^2)^2 + uy*(-(pixy*uy) + (ux*(piyy + piyy*ux^2 - 2*pixy*ux*uy + pizz*(1 + ux^2 + uy^2)))/(1 + uy^2)))/(1 + ux^2 + uy^2)^(3/2),
        (2*((piyy + pizz*ux^2 + piyy*ux^2)*uy + pixy*(ux - ux*uy^2)))/(1 + uy^2)^2,
        (piB + dtp*T)*ux,
        (2*ux*uy*(-(piyy*tauS*(1 + ux^2)) + pixy*tauS*ux*uy + (1 + ux^2 + uy^2)*visc))/(1 + ux^2 + uy^2),
        0,
        (-(piyy*tauS*ux^2*(1 + ux^2)*(1 + 2*uy^2)) + pixy*tauS*ux*uy*(-1 + (-1 + 2*ux^2)*uy^2) - (1 + ux^2 + uy^2)*(pizz*tauS*ux^2*uy^2 - (1 + ux^2)*(1 + uy^2)*visc))/((1 + uy^2)*(1 + ux^2 + uy^2)),
        0,
        -((ux*(1 + ux^2))/((1 + uy^2)*sqrt(1 + ux^2 + uy^2))),
        -((1 + ux^2)/(1 + uy^2)),
        0,
        tauS*ux,
        0,
        0,
        0,
        -((ux*sqrt(1 + ux^2 + uy^2))/(1 + uy^2)),
        -((1 + ux^2 + uy^2)/(1 + uy^2)),
        0,
        0,
        (tauS*ux)/tau^2,
        0,
        0,
        (uy*(1 + (2*ux^2)/(1 + uy^2)))/sqrt(1 + ux^2 + uy^2),
        (2*ux*uy)/(1 + uy^2),
        1,
        0,
        0,
        tauS*ux,
        0,
        ux*sqrt(1 + ux^2 + uy^2),
        1 + ux^2,
        ux*uy,
        0,
        0,
        0,
        tauB*ux))

    Ay = SMatrix{7,7}((
        (dtp + dtdtp*T)*uy*sqrt(1 + ux^2 + uy^2),
        (dtp + dtdtp*T)*ux*uy,
        dtp + dtp*uy^2 + dtdtp*T*uy^2,
        0,
        0,
        0,
        0,
        (pixy*(1 + uy^2) + ux*uy*(-piyy + piB*(1 + ux^2 + uy^2) + dtp*T*(1 + ux^2 + uy^2)))/(1 + ux^2 + uy^2)^(3/2),
        (piB + dtp*T)*uy,
        0,
        (-2*tauS*uy^2*(pixy - piyy*ux*uy + pixy*uy^2))/(1 + ux^2 + uy^2),
        0,
        (-(pixy*tauS*ux*uy*(1 + 2*uy^2)) + visc + uy^4*(pizz*tauS + visc) + uy^2*(pizz*tauS + piyy*tauS + 2*visc) + ux^2*(visc + uy^2*(pizz*tauS + 2*piyy*tauS + visc)))/(1 + ux^2 + uy^2),
        0,
        (-(uy*(pixy*ux + piyy*uy)) + piyy*(1 + ux^2 + uy^2) + (piB + dtp*T)*uy^2*(1 + ux^2 + uy^2) + (piB + dtp*T)*(1 + ux^2 + uy^2)^2)/(1 + ux^2 + uy^2)^(3/2),
        (piB + dtp*T)*ux,
        2*(piB + dtp*T)*uy,
        (-6*piyy*tauS*(1 + ux^2)*uy^2 + 6*pixy*tauS*ux*uy^3 + 4*(1 + uy^2)*(1 + ux^2 + uy^2)*visc)/(3*(1 + ux^2 + uy^2)),
        (-2*visc)/(3*tau^2),
        -1/3*(uy*(3*piyy*tauS*ux*(1 + ux^2)*(1 + 2*uy^2) + 3*pixy*tauS*uy*(1 + (1 - 2*ux^2)*uy^2) + ux*(1 + ux^2 + uy^2)*(3*pizz*tauS*uy^2 - (1 + uy^2)*visc)))/((1 + uy^2)*(1 + ux^2 + uy^2)),
        zeta,
        uy/sqrt(1 + ux^2 + uy^2),
        0,
        1,
        tauS*uy,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        (tauS*uy)/tau^2,
        0,
        0,
        ux/sqrt(1 + ux^2 + uy^2),
        1,
        0,
        0,
        0,
        tauS*uy,
        0,
        uy*sqrt(1 + ux^2 + uy^2),
        ux*uy,
        1 + uy^2,
        0,
        0,
        0,
        tauB*uy))

    source = SVector{7}((
        (dtp*T - piyy*ux^2 + dtp*T*ux^2 + 2*pixy*ux*uy + piyy*uy^2 + 2*dtp*T*uy^2 + dtp*T*ux^2*uy^2 + dtp*T*uy^4 + pizz*(1 - ux^2 + uy^2) + piB*(1 + uy^2)*(1 + ux^2 + uy^2))/(tau*(1 + uy^2)),
        (-(piyy*ux) + dtp*T*ux - piyy*ux^3 + dtp*T*ux^3 + pixy*uy + 2*pixy*ux^2*uy + 2*dtp*T*ux*uy^2 + dtp*T*ux^3*uy^2 + pixy*uy^3 + dtp*T*ux*uy^4 - pizz*ux*(1 + ux^2 + uy^2) + piB*ux*(1 + uy^2)*(1 + ux^2 + uy^2))/(tau*(1 + uy^2)*sqrt(1 + ux^2 + uy^2)),
        (pixy*ux + uy*(piyy + piB*(1 + ux^2 + uy^2) + dtp*T*(1 + ux^2 + uy^2)))/(tau*sqrt(1 + ux^2 + uy^2)),
        piyy - (2*(1 + uy^2)*sqrt(1 + ux^2 + uy^2)*visc)/(3*tau),
        (3*pizz*tau + 4*sqrt(1 + ux^2 + uy^2)*visc)/(3*tau^3),
        pixy - (2*ux*uy*sqrt(1 + ux^2 + uy^2)*visc)/(3*tau),
        piB + (sqrt(1 + ux^2 + uy^2)*zeta)/tau))

    iszero(dpp) && return (At, Ax, Ay, source)

    # the coefficient of dpp of  dpp pi^{mn} theta
    dAt = SMatrix{7,7}((
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        (piyy*ux)/sqrt(1 + ux^2 + uy^2),
        (pizz*ux)/(tau^2*sqrt(1 + ux^2 + uy^2)),
        (pixy*ux)/sqrt(1 + ux^2 + uy^2),
        0,
        0,
        0,
        0,
        (piyy*uy)/sqrt(1 + ux^2 + uy^2),
        (pizz*uy)/(tau^2*sqrt(1 + ux^2 + uy^2)),
        (pixy*uy)/sqrt(1 + ux^2 + uy^2),
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0))

    dAx = SMatrix{7,7}((
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        piyy,
        pizz/tau^2,
        pixy,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0))

    dAy = SMatrix{7,7}((
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        piyy,
        pizz/tau^2,
        pixy,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0))

    dsource = SVector{7}((
        0,
        0,
        0,
        (piyy*sqrt(1 + ux^2 + uy^2))/tau,
        (pizz*sqrt(1 + ux^2 + uy^2))/tau^3,
        (pixy*sqrt(1 + ux^2 + uy^2))/tau,
        0))

    return (At + dpp*dAt, Ax + dpp*dAx, Ay + dpp*dAy, source + dpp*dsource)
end
# end of generated code
