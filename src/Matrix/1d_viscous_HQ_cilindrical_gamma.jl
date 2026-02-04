# the convention here are T, ur,  \[Pi]phiphi, \[Pi]etaeta, \[Pi]B, mu, nur

#THIS IS THE MATRIX THAT DID NOT CREATE PROBLEMS WITH THE BUMP (GUBSER)

function matrxi1d_visc_HQ!(A_i,Source,ϕ,t,X,params)

    dpt = pressure_derivative(ϕ[1],Val(1),params.eos) #entropy
    dptt = pressure_derivative(ϕ[1],Val(2),params.eos)
        
    etaVisc=viscosity(ϕ[1],dpt,params.shear)
    tauS=τ_shear(ϕ[1],dpt,params.shear)
    tauB=τ_bulk(ϕ[1],dpt,dptt,params.bulk)
    zeta=bulk_viscosity(ϕ[1],dpt,params.bulk)
    
    thermo = thermodynamic(ϕ[1],0.0,params.eos.hadron_list)
    n=thermo.pressure
    dtn, dmn = thermo.pressure_derivative
    dmn+=0.0001
    #dttn, dtdmn, dmmn = thermodynamic(ϕ[1],ϕ[6],HadronResonaceGasNew()).pressure_hessian.* fmGeV^3
    # @show n, dtn, dmn

    # kappa = diffusion(ϕ[1],n,params.diffusion)
    # #@show ϕ[1]
    # #@show t, X[1], ϕ[1] #X[1] is the radius: here we print the time, the radius, and the temperature
    # tauDiff=τ_diffusion(ϕ[1],params.diffusion)
    kappa = diffusion_hadron(ϕ[1],0.0,params.eos,params.diffusion) #diffusion coefficient for hadrons
    tauDiff=τ_diffusion_hadron(ϕ[1],0.0,params.eos,params.diffusion) #tau diffusion for hadrons
    
    dmp = 0 #for now we don t have chemical potential
    dtdmp = 0 
    dmdmp = 0

    #actually our equations don t depend on p: we can just put as entry dpt instead, in any case it will not be used (but in the future maybe it will be )
    #(At,Ax, source)=one_d_viscous_HQ_matrix(ϕ,t,X[1],dpt,dpt,dptt,zeta,etaVisc,tauS,tauB,n,dtn,dmn,tauDiff,Ds)
    (At,Ax, source)=one_d_viscous_matrix(ϕ,t,X[1],dpt,dpt,dptt,dmp,dtdmp,dmdmp,zeta,etaVisc,tauS,tauB,n,dtn,dmn,tauDiff,kappa)

        
    Ainv= inv(At)


    A_mul_B!(A_i[1], Ainv,Ax)
    
    
    jgemvavx!(Source, Ainv,source)

        
    
    end 



function one_d_viscous_matrix(u,tau,R,p,dtp,dtdtp,dmp,dtdmp,dmdmp,zeta,visc,tauS,tauB,n,dtn,dmn,tauDiff,kappa)
    #these are the matrices from 1d viscous hydro
    At=SMatrix{7,7}(
        dtdtp*u[1] + (dtp + dtdtp*u[1])*^(u[2],2)

    ,(dtp + dtdtp*u[1])*u[2]*sqrt(1 + ^(u[2],2))

    ,0

    ,0

    ,0

    ,0

    ,0

    ,2*u[2]*(dtp*u[1] + u[5] - u[3] - u[4])

    ,((1 + 2*^(u[2],2))*(dtp*u[1] + u[5] - u[3] - u[4]))/sqrt(1 + ^(u[2],2))

    ,(-2*u[2]*(visc - 2*tauS*u[3]))/(3. *^(R,2)*sqrt(1 + ^(u[2],2)))

    ,(-2*u[2]*(visc - 2*tauS*u[4]))/(3. *^(tau,2)*sqrt(1 + ^(u[2],2)))

    ,(u[2]*zeta)/sqrt(1 + ^(u[2],2))

    ,0

    ,0

    ,-^(u[2],2)

    ,-(u[2]*sqrt(1 + ^(u[2],2)))

    ,(tauS*sqrt(1 + ^(u[2],2)))/^(R,2)

    ,0

    ,0

    ,0

    ,0

    ,-^(u[2],2)

    ,-(u[2]*sqrt(1 + ^(u[2],2)))

    ,0

    ,(tauS*sqrt(1 + ^(u[2],2)))/^(tau,2)

    ,0

    ,0

    ,0

    ,^(u[2],2)

    ,u[2]*sqrt(1 + ^(u[2],2))

    ,0

    ,0

    ,tauB*sqrt(1 + ^(u[2],2))

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,n*sqrt(1 + ^(u[2],2))

    ,kappa*u[2]*sqrt(1 + ^(u[2],2))

    ,0

    ,0

    ,0

    ,0

    ,0

    ,u[2]/sqrt(1 + ^(u[2],2))

    ,tauDiff*sqrt(1 + ^(u[2],2))
    )
        #########################################################################################################################################################################
        
    Ax=SMatrix{7,7}(
        (dtp + dtdtp*u[1])*u[2]*sqrt(1 + ^(u[2],2))

        ,dtp + (dtp + dtdtp*u[1])*^(u[2],2)

        ,0

        ,0

        ,0

        ,(-3*dtp*tauB*tauS*(dtdtp*u[1]*(u[7] + u[6]*n*u[2]*(1 + ^(u[2],2))) + dtn*u[6]*u[2]*(1 + ^(u[2],2))*(-(dtp*u[1]) - u[5] + u[3] + u[4])))/((1 + ^(u[2],2))*(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))

        ,(3*dtdtp*dtp*u[7]*u[1]*tauB*tauDiff*tauS*u[2])/(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4])))))

        ,((1 + 2*^(u[2],2))*(dtp*u[1] + u[5] - u[3] - u[4]))/sqrt(1 + ^(u[2],2))

        ,2*u[2]*(dtp*u[1] + u[5] - u[3] - u[4])

        ,(-2*(visc - 2*tauS*u[3]))/(3. *^(R,2))

        ,(-2*(visc - 2*tauS*u[4]))/(3. *^(tau,2))

        ,zeta

        ,(3*dtp*u[7]*tauB*tauS*u[2]*(dtp*u[1] + u[5] - u[3] - u[4]) - 3*dtn*u[6]*tauB*tauS*(1 + ^(u[2],2))*^(-(dtp*u[1]) - u[5] + u[3] + u[4],2) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(-(u[7]*u[2]) + u[6]*n*(1 + ^(u[2],2))) + 3*u[6]*n*tauB*tauS*(1 + ^(u[2],2))*(u[5] - u[3] - u[4]) + u[7]*u[2]*(4*tauB*visc + 3*tauS*zeta - tauB*tauS*(3*u[5] + u[3] + u[4]))))/((1 + ^(u[2],2))*(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))

        ,(u[7]*tauDiff*^(u[2],2)*(-3*^(dtp,2)*u[1]*tauB*tauS + 3*dtp*tauB*tauS*(-u[5] + u[3] + u[4]) + dtdtp*u[1]*(-4*tauB*visc - 3*tauS*zeta + 4*tauB*tauS*(u[3] + u[4]))))/((1 + ^(u[2],2))*(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))

        ,-(u[2]*sqrt(1 + ^(u[2],2)))

        ,-1 - ^(u[2],2)

        ,(tauS*u[2])/^(R,2)

        ,0

        ,0

        ,(3*tauB*tauS*(dtdtp*u[1]*(u[7] + u[6]*n*u[2]*(1 + ^(u[2],2))) + dtn*u[6]*u[2]*(1 + ^(u[2],2))*(-(dtp*u[1]) - u[5] + u[3] + u[4])))/((1 + ^(u[2],2))*(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))

        ,(-3*dtdtp*u[7]*u[1]*tauB*tauDiff*tauS*u[2])/(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4])))))

        ,-(u[2]*sqrt(1 + ^(u[2],2)))

        ,-1 - ^(u[2],2)

        ,0

        ,(tauS*u[2])/^(tau,2)

        ,0

        ,(3*tauB*tauS*(dtdtp*u[1]*(u[7] + u[6]*n*u[2]*(1 + ^(u[2],2))) + dtn*u[6]*u[2]*(1 + ^(u[2],2))*(-(dtp*u[1]) - u[5] + u[3] + u[4])))/((1 + ^(u[2],2))*(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))

        ,(-3*dtdtp*u[7]*u[1]*tauB*tauDiff*tauS*u[2])/(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4])))))

        ,u[2]*sqrt(1 + ^(u[2],2))

        ,1 + ^(u[2],2)

        ,0

        ,0

        ,tauB*u[2]

        ,(-3*tauB*tauS*(dtdtp*u[1]*(u[7] + u[6]*n*u[2]*(1 + ^(u[2],2))) + dtn*u[6]*u[2]*(1 + ^(u[2],2))*(-(dtp*u[1]) - u[5] + u[3] + u[4])))/((1 + ^(u[2],2))*(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))

        ,(3*dtdtp*u[7]*u[1]*tauB*tauDiff*tauS*u[2])/(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4])))))

        ,0

        ,0

        ,0

        ,0

        ,0

        ,n*u[2]

        ,kappa*(1 + ^(u[2],2))

        ,0

        ,0

        ,0

        ,0

        ,0

        ,1

        ,tauDiff*u[2]
    )
        #########################################################################################################################################################################
    
    source=SVector{7}(
        (dtp*u[1]*(R + R*^(u[2],2) + tau*u[2]*sqrt(1 + ^(u[2],2))) + tau*u[2]*sqrt(1 + ^(u[2],2))*(u[5] - u[3] - u[4]) + R*((1 + ^(u[2],2))*u[5] + u[4] - ^(u[2],2)*(u[3] + u[4])))/(R*tau)

        ,(dtp*u[1]*u[2]*(tau*u[2] + R*sqrt(1 + ^(u[2],2))) + R*u[2]*sqrt(1 + ^(u[2],2))*(u[5] - u[3] - u[4]) + tau*(-2*u[3] + ^(u[2],2)*(u[5] - u[3] - u[4]) - u[4]))/(R*tau)

        ,(3*R*u[3] - (2*R*sqrt(1 + ^(u[2],2))*(visc - 2*tauS*u[3]))/tau + 4*u[2]*(visc + tauS*u[3]))/(3. *^(R,3))

        ,(3*tau*u[4] - (2*tau*u[2]*(visc - 2*tauS*u[4]))/R + 4*sqrt(1 + ^(u[2],2))*(visc + tauS*u[4]))/(3. *^(tau,3))

        ,(u[2]*zeta)/R + (sqrt(1 + ^(u[2],2))*zeta)/tau + u[5]

        ,(3*dtp*tauB*tauS*u[2]*(dtp*u[7]*u[1]*(tau*^(u[2],3)*sqrt(1 + ^(u[2],2)) + R*(-1 + ^(u[2],4))) - u[6]*n*u[2]*(1 + ^(u[2],2))*(tau*u[2]*sqrt(1 + ^(u[2],2))*(2*u[3] + u[4]) + R*(1 + ^(u[2],2))*(u[3] + 2*u[4])) + u[7]*(tau*u[2]*sqrt(1 + ^(u[2],2))*(-2*u[3] + ^(u[2],2)*(u[5] - u[3] - u[4]) - u[4]) + R*(1 + ^(u[2],2))*((-1 + ^(u[2],2))*u[5] - u[4] - ^(u[2],2)*(u[3] + u[4])))) + dtn*u[6]*(1 + ^(u[2],2))*(3*^(dtp,2)*^(u[1],2)*tauB*tauS*(1 + ^(u[2],2))*(R + R*^(u[2],2) + tau*u[2]*sqrt(1 + ^(u[2],2))) - tau*u[2]*sqrt(1 + ^(u[2],2))*(3*tauS*^(u[2],2)*zeta*(2*u[3] + u[4]) - 3*tauB*tauS*(u[5] - u[3] - u[4])*(u[5] + 3*u[3] + u[4]) + tauB*^(u[2],2)*(2*visc*(3*u[5] + u[3] - u[4]) - tauS*(3*^(u[5],2) + 6*u[5]*u[3] - ^(u[3],2) + ^(u[4],2)))) - 3*dtp*u[1]*(2*tau*tauB*u[2]*sqrt(1 + ^(u[2],2))*(-(tauS*(u[5] + u[3])) + ^(u[2],2)*(visc - tauS*(u[5] + u[3]))) + R*(-(tau*tauS*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[5]) + tauB*(-2*tauS*u[5] + tauS*u[3] + 2*^(u[2],4)*(visc - tauS*(u[5] + u[4])) + ^(u[2],2)*(2*visc + tauS*(-4*u[5] + u[3] - 2*u[4]) + tau*sqrt(1 + ^(u[2],2))*(u[3] + u[4]))))) - R*(3*tauS*^(u[2],2)*(tau*sqrt(1 + ^(u[2],2))*u[5]*(-u[5] + u[3] + u[4]) + (1 + ^(u[2],2))*zeta*(u[3] + 2*u[4])) + tauB*(-3*tauS*(u[5] - u[3] - u[4])*(u[5] + u[4]) + ^(u[2],4)*(2*visc*(3*u[5] - u[3] + u[4]) - tauS*(3*^(u[5],2) + ^(u[3],2) + 6*u[5]*u[4] - ^(u[4],2))) + ^(u[2],2)*(2*visc*(3*u[5] - u[3] + u[4]) + 3*tau*sqrt(1 + ^(u[2],2))*(u[5] - u[3] - u[4])*(u[3] + u[4]) - tauS*(6*^(u[5],2) - 3*u[5]*(u[3] - 2*u[4]) + (u[3] - 4*u[4])*(u[3] + u[4])))))) - dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*^(1 + ^(u[2],2),2)*(u[7]*(R*u[2] + tau*sqrt(1 + ^(u[2],2))) + u[6]*n*(R + R*^(u[2],2) + tau*u[2]*sqrt(1 + ^(u[2],2)))) - 3*u[6]*n*(1 + ^(u[2],2))*(tau*tauB*u[2]*sqrt(1 + ^(u[2],2))*(-(tauS*(u[5] + u[3])) + ^(u[2],2)*(2*visc - tauS*(u[5] + u[3]))) + R*(-(tau*tauS*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[5]) + tauB*(tauS*(-u[5] + u[3] + u[4]) + ^(u[2],4)*(2*visc - tauS*(u[5] + u[4])) + ^(u[2],2)*(2*visc + tauS*(-2*u[5] + u[3]) + tau*sqrt(1 + ^(u[2],2))*(u[3] + u[4]))))) + u[7]*(R*u[2]*(-2*tauB*(1 + 3*^(u[2],2) + 2*^(u[2],4))*visc + 3*tauS*(zeta - ^(u[2],4)*zeta + tau*sqrt(1 + ^(u[2],2))*u[5]) - 3*tau*tauB*sqrt(1 + ^(u[2],2))*(u[3] + u[4]) + tauB*tauS*(1 + ^(u[2],2))*(3*(1 + ^(u[2],2))*u[5] - 4*u[3] - u[4] + ^(u[2],2)*(u[3] + u[4]))) + tau*sqrt(1 + ^(u[2],2))*(-3*tauS*^(u[2],4)*zeta + tauB*(3*tauS*(u[5] + u[3]) - 3*^(u[2],2)*(2*visc - 2*tauS*u[5] + tauS*u[4]) + ^(u[2],4)*(-4*visc + tauS*(3*u[5] + u[3] + u[4])))))))/(R*tau*^(1 + ^(u[2],2),1.5)*(3*dtp*tauB*tauS*^(u[2],2)*(dtp*u[1] + u[5] - u[3] - u[4]) - dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))

        ,(u[7]*(-3*dtp*tauB*tauS*^(u[2],2)*(dtp*u[1]*(tau*tauDiff*u[2] + R*(tau + tauDiff*sqrt(1 + ^(u[2],2)))) + tau*tauDiff*u[2]*(u[5] + u[3]) + R*tau*(u[5] - u[3] - u[4]) + R*tauDiff*sqrt(1 + ^(u[2],2))*(u[5] + u[4])) + dtdtp*u[1]*(3*dtp*R*u[1]*tau*tauB*tauS*(1 + ^(u[2],2)) + tau*tauDiff*u[2]*(-3*tauS*^(u[2],2)*zeta - 3*tauB*tauS*(2*u[3] + u[4]) + tauB*^(u[2],2)*(2*visc + tauS*(-2*u[3] + u[4]))) + R*(-3*tau*tauS*^(u[2],2)*(zeta + tauDiff*u[5]) + tauDiff*^(u[2],2)*sqrt(1 + ^(u[2],2))*(2*tauB*visc - 3*tauS*zeta + tauB*tauS*(u[3] - 2*u[4])) + tau*tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + 3*tauDiff*(u[3] + u[4]) + tauS*(3*u[5] + u[3] + u[4])))))))/(R*tau*(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))
    )
        
    return (At,Ax, source)
end

function A_mul_B!(C, A, B)
    @turbo for n ∈ indices((C,B), 2), m ∈ indices((C,A), 1)
        Cmn = zero(eltype(C))
        for k ∈ indices((A,B), (2,1))
            Cmn += A[m,k] * B[k,n]
        end
        C[m,n] = Cmn
    end
end
