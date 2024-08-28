# the convention here are T, ur,  \[Pi]phiphi, \[Pi]etaeta, \[Pi]B, mu, nur

function matrxi1d_visc_HQ!(A_i,Source,ϕ,t,X,params)
    
    #@show t, ϕ
    dpt = pressure_derivative(ϕ[1],Val(1),params.eos) #entropy

    

    dptt = pressure_derivative(ϕ[1],Val(2),params.eos)
        
    if params.shear.ηs == 0 #zero bulk case
        tauS=1
        etaVisc=0
    
    else #non 
    etaVisc=viscosity(ϕ[1],dpt,params.shear)
    
    tauS=τ_shear(ϕ[1],dpt,params.shear)
    end
    
    if params.bulk.ζs == 0 #zero bulk case
        tauB=1
        zeta=0
    
    else #non zero bulk case    
    tauB=τ_bulk(ϕ[1],dpt,dptt,params.bulk)
    zeta=bulk_viscosity(ϕ[1],dpt,params.bulk)
    end
    #end
    
    #n,dtn,dmn,dtmn = federica(ϕ[1],ϕ[6],params.eos)
    thermo = thermodynamic(ϕ[1],ϕ[6],params.eos.hadron_list)
    n=thermo.pressure
    dtn, dmn = thermo.pressure_derivative
    dmn+=0.0001
    #dtn+=0.0001
    #dttn, dtdmn, dmmn = thermodynamic(ϕ[1],ϕ[6],HadronResonaceGasNew()).pressure_hessian.* fmGeV^3
    
    Ds = diffusion(ϕ[1],n,params.diffusion)
    #@show ϕ[1]
    #@show t, X[1], ϕ[1] #X[1] is the radius: here we print the time, the radius, and the temperature
    tauDiff=τ_diffusion(ϕ[1],params.diffusion)
    #@show t
    dmp = 0 #for now we don t have chemical potential in the eos
    dtdmp = 0 
    dmdmp = 0

    #actually our equations don t depend on p: we can just put as entry dpt instead, in any case it will not be used (but in the future maybe it will be )

    #(At,Ax, source)=one_d_viscous_matrix(ϕ,t,X[1],dpt,dpt,dptt,dmp,dtdmp,dmdmp,zeta,etaVisc,tauS,tauB,n,dtn,dmn,tauDiff,Ds)
    (At,Ax, source)=one_d_viscous_matrix(ϕ,t,X[1],dpt,dpt,dptt,zeta,etaVisc,tauS,tauB,n,dtn,dmn,tauDiff,Ds)
        
    Ainv= inv(At)


    A_mul_B!(A_i[1], Ainv,Ax)
    
    
    jgemvavx!(Source, Ainv,source)

        
    
    end 

    
    function one_d_viscous_matrix(u,tau,R,p,dtp,dtdtp,zeta,visc,tauS,tauB,n,dtn,dmn,tauDiff,Ds)
    
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
            
            ,(-2*u[2]*visc)/(3*^(R,2)*sqrt(1 + ^(u[2],2)))
            
            ,(-2*u[2]*visc)/(3*^(tau,2)*sqrt(1 + ^(u[2],2)))
            
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
            
            ,dmn*sqrt(1 + ^(u[2],2))
            
            ,Ds*n*u[2]*sqrt(1 + ^(u[2],2))
            
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
            
            ,(-3*dtp*tauB*tauS*(dtdtp*u[1]*(u[7] + n*u[2]*(1 + ^(u[2],2))) + dtn*(u[2] + ^(u[2],3))*(-(dtp*u[1]) - u[5] + u[3] + u[4])))/((1 + ^(u[2],2))*(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta - tauB*(3*tauS*(-u[5] + u[3] + u[4]) + ^(u[2],2)*(4*visc + 3*tauS*(-u[5] + u[3] + u[4]))))))
            
            ,(3*dtdtp*dtp*u[7]*u[1]*tauB*tauDiff*tauS*u[2])/(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta - tauB*(3*tauS*(-u[5] + u[3] + u[4]) + ^(u[2],2)*(4*visc + 3*tauS*(-u[5] + u[3] + u[4])))))
            
            ,((1 + 2*^(u[2],2))*(dtp*u[1] + u[5] - u[3] - u[4]))/sqrt(1 + ^(u[2],2))
            
            ,2*u[2]*(dtp*u[1] + u[5] - u[3] - u[4])
            
            ,(-2*visc)/(3*^(R,2))
            
            ,(-2*visc)/(3*^(tau,2))
            
            ,zeta
            
            ,(3*dtp*u[7]*tauB*tauS*u[2]*(dtp*u[1] + u[5] - u[3] - u[4]) - 3*dtn*tauB*tauS*(1 + ^(u[2],2))*^(-(dtp*u[1]) - u[5] + u[3] + u[4],2) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(n - u[7]*u[2] + n*^(u[2],2)) + 3*n*tauB*tauS*(1 + ^(u[2],2))*(u[5] - u[3] - u[4]) + u[7]*u[2]*(4*tauB*visc + 3*tauS*zeta + 3*tauB*tauS*(-u[5] + u[3] + u[4]))))/((1 + ^(u[2],2))*(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta - tauB*(3*tauS*(-u[5] + u[3] + u[4]) + ^(u[2],2)*(4*visc + 3*tauS*(-u[5] + u[3] + u[4]))))))
            
            ,-((u[7]*tauDiff*^(u[2],2)*(3*^(dtp,2)*u[1]*tauB*tauS + dtdtp*u[1]*(4*tauB*visc + 3*tauS*zeta) + 3*dtp*tauB*tauS*(u[5] - u[3] - u[4])))/((1 + ^(u[2],2))*(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta - tauB*(3*tauS*(-u[5] + u[3] + u[4]) + ^(u[2],2)*(4*visc + 3*tauS*(-u[5] + u[3] + u[4])))))))
            
            ,-(u[2]*sqrt(1 + ^(u[2],2)))
            
            ,-1 - ^(u[2],2)
            
            ,(tauS*u[2])/^(R,2)
            
            ,0
            
            ,0
            
            ,(3*tauB*tauS*(dtdtp*u[1]*(u[7] + n*u[2]*(1 + ^(u[2],2))) + dtn*(u[2] + ^(u[2],3))*(-(dtp*u[1]) - u[5] + u[3] + u[4])))/((1 + ^(u[2],2))*(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta - tauB*(3*tauS*(-u[5] + u[3] + u[4]) + ^(u[2],2)*(4*visc + 3*tauS*(-u[5] + u[3] + u[4]))))))
            
            ,(-3*dtdtp*u[7]*u[1]*tauB*tauDiff*tauS*u[2])/(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta - tauB*(3*tauS*(-u[5] + u[3] + u[4]) + ^(u[2],2)*(4*visc + 3*tauS*(-u[5] + u[3] + u[4])))))
            
            ,-(u[2]*sqrt(1 + ^(u[2],2)))
            
            ,-1 - ^(u[2],2)
            
            ,0
            
            ,(tauS*u[2])/^(tau,2)
            
            ,0
            
            ,(3*tauB*tauS*(dtdtp*u[1]*(u[7] + n*u[2]*(1 + ^(u[2],2))) + dtn*(u[2] + ^(u[2],3))*(-(dtp*u[1]) - u[5] + u[3] + u[4])))/((1 + ^(u[2],2))*(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta - tauB*(3*tauS*(-u[5] + u[3] + u[4]) + ^(u[2],2)*(4*visc + 3*tauS*(-u[5] + u[3] + u[4]))))))
            
            ,(-3*dtdtp*u[7]*u[1]*tauB*tauDiff*tauS*u[2])/(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta - tauB*(3*tauS*(-u[5] + u[3] + u[4]) + ^(u[2],2)*(4*visc + 3*tauS*(-u[5] + u[3] + u[4])))))
            
            ,u[2]*sqrt(1 + ^(u[2],2))
            
            ,1 + ^(u[2],2)
            
            ,0
            
            ,0
            
            ,tauB*u[2]
            
            ,(-3*tauB*tauS*(dtdtp*u[1]*(u[7] + n*u[2]*(1 + ^(u[2],2))) + dtn*(u[2] + ^(u[2],3))*(-(dtp*u[1]) - u[5] + u[3] + u[4])))/((1 + ^(u[2],2))*(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta - tauB*(3*tauS*(-u[5] + u[3] + u[4]) + ^(u[2],2)*(4*visc + 3*tauS*(-u[5] + u[3] + u[4]))))))
            
            ,(3*dtdtp*u[7]*u[1]*tauB*tauDiff*tauS*u[2])/(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta - tauB*(3*tauS*(-u[5] + u[3] + u[4]) + ^(u[2],2)*(4*visc + 3*tauS*(-u[5] + u[3] + u[4])))))
            
            ,0
            
            ,0
            
            ,0
            
            ,0
            
            ,0
            
            ,dmn*u[2]
            
            ,Ds*n*(1 + ^(u[2],2))
            
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
            
            ,(4*tau*u[2]*visc - 2*R*sqrt(1 + ^(u[2],2))*visc + 3*R*tau*u[3])/(3*^(R,3)*tau)
            
            ,(-2*tau*u[2]*visc + 4*R*sqrt(1 + ^(u[2],2))*visc + 3*R*tau*u[4])/(3*R*^(tau,3))
            
            ,(u[2]*zeta)/R + (sqrt(1 + ^(u[2],2))*zeta)/tau + u[5]
            
            ,(3*dtp*tauB*tauS*u[2]*(dtp*u[7]*u[1]*(tau*^(u[2],3)*sqrt(1 + ^(u[2],2)) + R*(-1 + ^(u[2],4))) - n*u[2]*(1 + ^(u[2],2))*(tau*u[2]*sqrt(1 + ^(u[2],2))*(2*u[3] + u[4]) + R*(1 + ^(u[2],2))*(u[3] + 2*u[4])) + u[7]*(tau*u[2]*sqrt(1 + ^(u[2],2))*(-2*u[3] + ^(u[2],2)*(u[5] - u[3] - u[4]) - u[4]) + R*(1 + ^(u[2],2))*((-1 + ^(u[2],2))*u[5] - u[4] - ^(u[2],2)*(u[3] + u[4])))) - dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*^(1 + ^(u[2],2),2)*(u[7]*R*u[2] + u[7]*tau*sqrt(1 + ^(u[2],2)) + n*tau*u[2]*sqrt(1 + ^(u[2],2)) + n*R*(1 + ^(u[2],2))) - 3*n*(1 + ^(u[2],2))*(tau*tauB*u[2]*sqrt(1 + ^(u[2],2))*(-(tauS*(u[5] + u[3])) + ^(u[2],2)*(2*visc - tauS*(u[5] + u[3]))) + R*(-(tau*tauS*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[5]) + tauB*(tauS*(-u[5] + u[3] + u[4]) + ^(u[2],4)*(2*visc - tauS*(u[5] + u[4])) + ^(u[2],2)*(2*visc + tauS*(-2*u[5] + u[3]) + tau*sqrt(1 + ^(u[2],2))*(u[3] + u[4]))))) - u[7]*(R*u[2]*(tauB*(2 + 6*^(u[2],2) + 4*^(u[2],4))*visc + 3*tauS*(-1 + ^(u[2],4))*zeta - 3*tau*tauS*sqrt(1 + ^(u[2],2))*u[5] + 3*tau*tauB*sqrt(1 + ^(u[2],2))*(u[3] + u[4]) - 3*tauB*tauS*(1 + ^(u[2],2))*((1 + ^(u[2],2))*u[5] + u[4] - ^(u[2],2)*(u[3] + u[4]))) + tau*sqrt(1 + ^(u[2],2))*(3*tauS*^(u[2],4)*zeta + tauB*(-3*tauS*(u[5] + u[3]) + 3*^(u[2],2)*(2*visc - 2*tauS*u[5] + tauS*u[4]) + ^(u[2],4)*(4*visc + 3*tauS*(-u[5] + u[3] + u[4])))))) + dtn*(1 + ^(u[2],2))*(3*^(dtp,2)*^(u[1],2)*tauB*tauS*(1 + ^(u[2],2))*(R + R*^(u[2],2) + tau*u[2]*sqrt(1 + ^(u[2],2))) - tau*u[2]*sqrt(1 + ^(u[2],2))*(3*tauS*^(u[2],2)*zeta*(2*u[3] + u[4]) - 3*tauB*tauS*(u[5] - u[3] - u[4])*(u[5] + 3*u[3] + u[4]) + tauB*^(u[2],2)*(2*visc*(3*u[5] + u[3] - u[4]) - 3*tauS*(u[5] - u[3] - u[4])*(u[5] + 3*u[3] + u[4]))) - 3*dtp*u[1]*(2*tau*tauB*u[2]*sqrt(1 + ^(u[2],2))*(-(tauS*(u[5] + u[3])) + ^(u[2],2)*(visc - tauS*(u[5] + u[3]))) + R*(-(tau*tauS*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[5]) + tauB*(-2*tauS*u[5] + tauS*u[3] + 2*^(u[2],4)*(visc - tauS*(u[5] + u[4])) + ^(u[2],2)*(2*visc + tauS*(-4*u[5] + u[3] - 2*u[4]) + tau*sqrt(1 + ^(u[2],2))*(u[3] + u[4]))))) - R*(3*tauS*^(u[2],2)*(tau*sqrt(1 + ^(u[2],2))*u[5]*(-u[5] + u[3] + u[4]) + (1 + ^(u[2],2))*zeta*(u[3] + 2*u[4])) + tauB*(-3*tauS*(u[5] - u[3] - u[4])*(u[5] + u[4]) + ^(u[2],4)*(2*visc*(3*u[5] - u[3] + u[4]) - 3*tauS*(u[5] - u[3] - u[4])*(u[5] + u[3] + 3*u[4])) + ^(u[2],2)*(2*visc*(3*u[5] - u[3] + u[4]) - 3*(u[5] - u[3] - u[4])*(-(tau*sqrt(1 + ^(u[2],2))*(u[3] + u[4])) + tauS*(2*u[5] + u[3] + 4*u[4])))))))/(R*tau*^(1 + ^(u[2],2),1.5)*(3*dtp*tauB*tauS*^(u[2],2)*(dtp*u[1] + u[5] - u[3] - u[4]) + dtdtp*u[1]*(-3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) + 3*tauS*^(u[2],2)*zeta + 3*tauB*tauS*(-u[5] + u[3] + u[4]) + tauB*^(u[2],2)*(4*visc + 3*tauS*(-u[5] + u[3] + u[4])))))
            
            ,(u[7]*(3*dtp*tauB*tauS*^(u[2],2)*(dtp*u[1]*(tau*tauDiff*u[2] + R*(tau + tauDiff*sqrt(1 + ^(u[2],2)))) + tau*tauDiff*u[2]*(u[5] + u[3]) + R*tau*(u[5] - u[3] - u[4]) + R*tauDiff*sqrt(1 + ^(u[2],2))*(u[5] + u[4])) + dtdtp*u[1]*(-3*dtp*R*u[1]*tau*tauB*tauS*(1 + ^(u[2],2)) + tau*tauDiff*u[2]*(3*tauS*^(u[2],2)*zeta + 3*tauB*tauS*(2*u[3] + u[4]) + tauB*^(u[2],2)*(-2*visc + 3*tauS*(2*u[3] + u[4]))) + R*(tauDiff*^(u[2],2)*sqrt(1 + ^(u[2],2))*(-2*tauB*visc + 3*tauS*zeta + 3*tauB*tauS*(u[3] + 2*u[4])) + tau*(3*tauS*^(u[2],2)*(zeta + tauDiff*u[5]) + 3*tauB*tauS*(-u[5] + u[3] + u[4]) + tauB*^(u[2],2)*(4*visc - 3*(tauS*(u[5] - u[3] - u[4]) + tauDiff*(u[3] + u[4]))))))))/(R*tau*(3*dtp*tauB*tauS*^(u[2],2)*(dtp*u[1] + u[5] - u[3] - u[4]) + dtdtp*u[1]*(-3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) + 3*tauS*^(u[2],2)*zeta + 3*tauB*tauS*(-u[5] + u[3] + u[4]) + tauB*^(u[2],2)*(4*visc + 3*tauS*(-u[5] + u[3] + u[4])))))
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

    #=
#gubser
    function one_d_viscous_matrix(u,tau,R,p,dtp,dtdtp,dmp,dtdmp,dmdmp,zeta,visc,tauS,tauB,n,dtn,dmn,tauDiff,Ds)
        #these are the matrices from 1d viscous hydro
        At=SMatrix{7,7}(
            dtdmp*u[6] + dtdtp*u[1] + (dtp + dtdmp*u[6] + dtdtp*u[1])*^(u[2],2)
    
            ,(dtp + dtdmp*u[6] + dtdtp*u[1])*u[2]*sqrt(1 + ^(u[2],2))
    
            ,0
    
            ,0
    
            ,0
    
            ,0
    
            ,0
    
            ,2*u[2]*(dmp*u[6] + dtp*u[1] + u[5] - u[3] - u[4])
    
            ,((1 + 2*^(u[2],2))*(dmp*u[6] + dtp*u[1] + u[5] - u[3] - u[4]))/sqrt(1 + ^(u[2],2))
    
            ,(-2*u[2]*(visc - 2*tauS*u[3]))/(3*^(R,2)*sqrt(1 + ^(u[2],2)))
    
            ,(-2*u[2]*(visc - 2*tauS*u[4]))/(3*^(tau,2)*sqrt(1 + ^(u[2],2)))
    
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
    
            ,dmdmp*u[6] + dtdmp*u[1] + (dmp + dmdmp*u[6] + dtdmp*u[1])*^(u[2],2)
    
            ,(dmp + dmdmp*u[6] + dtdmp*u[1])*u[2]*sqrt(1 + ^(u[2],2))
    
            ,0
    
            ,0
    
            ,0
    
            ,-((3*dmp*dtdmp*dtn*u[6]*u[1]*tauB*tauS + 3*dtdmp*dtn*dtp*^(u[1],2)*tauB*tauS + 3*dmp*dtdmp*u[6]*u[7]*tauB*tauS*u[2] + 3*dmp*dtdtp*u[7]*u[1]*tauB*tauS*u[2] - 3*dtdmp*dtp*u[7]*u[1]*tauB*tauS*u[2] - 3*^(dmp,2)*dtn*u[6]*tauB*tauS*^(u[2],2) + 3*dmp*dtdmp*u[6]*n*tauB*tauS*^(u[2],2) - 3*dmp*dtn*dtp*u[1]*tauB*tauS*^(u[2],2) + 6*dmp*dtdmp*dtn*u[6]*u[1]*tauB*tauS*^(u[2],2) + 3*dmp*dtdtp*n*u[1]*tauB*tauS*^(u[2],2) - 3*dtdmp*dtp*n*u[1]*tauB*tauS*^(u[2],2) + 6*dtdmp*dtn*dtp*^(u[1],2)*tauB*tauS*^(u[2],2) - 3*^(dmp,2)*dtn*u[6]*tauB*tauS*^(u[2],4) + 3*dmp*dtdmp*u[6]*n*tauB*tauS*^(u[2],4) - 3*dmp*dtn*dtp*u[1]*tauB*tauS*^(u[2],4) + 3*dmp*dtdmp*dtn*u[6]*u[1]*tauB*tauS*^(u[2],4) + 3*dmp*dtdtp*n*u[1]*tauB*tauS*^(u[2],4) - 3*dtdmp*dtp*n*u[1]*tauB*tauS*^(u[2],4) + 3*dtdmp*dtn*dtp*^(u[1],2)*tauB*tauS*^(u[2],4) - 4*dtdmp*dtn*u[1]*tauB*^(u[2],2)*visc - 4*dtdmp*dtn*u[1]*tauB*^(u[2],4)*visc - 3*dtdmp*dtn*u[1]*tauS*^(u[2],2)*zeta - 3*dtdmp*dtn*u[1]*tauS*^(u[2],4)*zeta + 3*dtdmp*dtn*u[1]*tauB*tauS*u[5] - 3*dmp*dtn*tauB*tauS*^(u[2],2)*u[5] + 6*dtdmp*dtn*u[1]*tauB*tauS*^(u[2],2)*u[5] - 3*dmp*dtn*tauB*tauS*^(u[2],4)*u[5] + 3*dtdmp*dtn*u[1]*tauB*tauS*^(u[2],4)*u[5] - 3*dtdmp*dtn*u[1]*tauB*tauS*u[3] + 3*dmp*dtn*tauB*tauS*^(u[2],2)*u[3] - 2*dtdmp*dtn*u[1]*tauB*tauS*^(u[2],2)*u[3] + 3*dmp*dtn*tauB*tauS*^(u[2],4)*u[3] + dtdmp*dtn*u[1]*tauB*tauS*^(u[2],4)*u[3] + dtn*tauB*tauS*(1 + ^(u[2],2))*(3*dmp*^(u[2],2) + dtdmp*u[1]*(-3 + ^(u[2],2)))*u[4] - dmn*(1 + ^(u[2],2))*(3*dtdtp*dtp*^(u[1],2)*tauB*tauS - 3*^(dtp,2)*u[1]*tauB*tauS*^(u[2],2) + 3*dtdtp*dtp*^(u[1],2)*tauB*tauS*^(u[2],2) + 3*dmp*u[6]*tauB*tauS*(dtdmp*u[6] + dtdtp*u[1] + (-dtp + dtdmp*u[6] + dtdtp*u[1])*^(u[2],2)) - 4*dtdtp*u[1]*tauB*^(u[2],2)*visc - 3*dtdtp*u[1]*tauS*^(u[2],2)*zeta + 3*dtdtp*u[1]*tauB*tauS*u[5] - 3*dtp*tauB*tauS*^(u[2],2)*u[5] + 3*dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[5] - 3*dtdtp*u[1]*tauB*tauS*u[3] + 3*dtp*tauB*tauS*^(u[2],2)*u[3] + dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[3] + tauB*tauS*(3*dtp*^(u[2],2) + dtdtp*u[1]*(-3 + ^(u[2],2)))*u[4] + dtdmp*u[6]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))) + dmdmp*u[6]*(3*dmp*dtn*u[6]*tauB*tauS*^(1 + ^(u[2],2),2) - 3*dtp*tauB*tauS*u[2]*(u[7] + n*u[2]*(1 + ^(u[2],2))) + dtn*(1 + ^(u[2],2))*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))/(sqrt(1 + ^(u[2],2))*(3*dtdtp*dtp*^(u[1],2)*tauB*tauS - 3*^(dtp,2)*u[1]*tauB*tauS*^(u[2],2) + 3*dtdtp*dtp*^(u[1],2)*tauB*tauS*^(u[2],2) + 3*dmp*u[6]*tauB*tauS*(dtdmp*u[6] + dtdtp*u[1] + (-dtp + dtdmp*u[6] + dtdtp*u[1])*^(u[2],2)) - 4*dtdtp*u[1]*tauB*^(u[2],2)*visc - 3*dtdtp*u[1]*tauS*^(u[2],2)*zeta + 3*dtdtp*u[1]*tauB*tauS*u[5] - 3*dtp*tauB*tauS*^(u[2],2)*u[5] + 3*dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[5] - 3*dtdtp*u[1]*tauB*tauS*u[3] + 3*dtp*tauB*tauS*^(u[2],2)*u[3] + dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[3] + tauB*tauS*(3*dtp*^(u[2],2) + dtdtp*u[1]*(-3 + ^(u[2],2)))*u[4] + dtdmp*u[6]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4])))))))
    
            ,(u[2]*(-3*dmp*dtdmp*u[6]*u[7]*tauB*tauDiff*tauS*u[2] + 3*dmdmp*dtp*u[6]*u[7]*tauB*tauDiff*tauS*u[2] - 3*dmp*dtdtp*u[7]*u[1]*tauB*tauDiff*tauS*u[2] + 3*dtdmp*dtp*u[7]*u[1]*tauB*tauDiff*tauS*u[2] + Ds*n*sqrt(1 + ^(u[2],2))*(3*(dtdmp*u[6] + dtdtp*u[1])*tauB*tauS*(dmp*u[6] + dtp*u[1] + u[5] - u[3] - u[4]) + ^(u[2],2)*(-3*^(dtp,2)*u[1]*tauB*tauS + 3*dmp*u[6]*(-dtp + dtdmp*u[6] + dtdtp*u[1])*tauB*tauS + 3*dtp*tauB*tauS*(dtdmp*u[6]*u[1] + dtdtp*^(u[1],2) - u[5] + u[3] + u[4]) - (dtdmp*u[6] + dtdtp*u[1])*(4*tauB*visc + 3*tauS*zeta - tauB*tauS*(3*u[5] + u[3] + u[4]))))))/(3*(dtdmp*u[6] + dtdtp*u[1])*tauB*tauS*(dmp*u[6] + dtp*u[1] + u[5] - u[3] - u[4]) + ^(u[2],2)*(-3*^(dtp,2)*u[1]*tauB*tauS + 3*dmp*u[6]*(-dtp + dtdmp*u[6] + dtdtp*u[1])*tauB*tauS + 3*dtp*tauB*tauS*(dtdmp*u[6]*u[1] + dtdtp*^(u[1],2) - u[5] + u[3] + u[4]) - (dtdmp*u[6] + dtdtp*u[1])*(4*tauB*visc + 3*tauS*zeta - tauB*tauS*(3*u[5] + u[3] + u[4]))))
    
            ,0
    
            ,0
    
            ,0
    
            ,0
    
            ,0
    
            ,u[2]/sqrt(1 + ^(u[2],2))
    
            ,-tauDiff
        )
            #########################################################################################################################################################################
            
        Ax=SMatrix{7,7}(
            (dtp + dtdmp*u[6] + dtdtp*u[1])*u[2]*sqrt(1 + ^(u[2],2))
    
            ,dtp + (dtp + dtdmp*u[6] + dtdtp*u[1])*^(u[2],2)
            
            ,0
            
            ,0
            
            ,0
            
            ,(3*dtp*tauB*tauS*(-(dtdtp*u[1]*(u[7] + n*u[2] + n*^(u[2],3))) - dtdmp*u[6]*(u[7] + n*u[2]*(1 + ^(u[2],2))) + dtn*(u[2] + ^(u[2],3))*(dmp*u[6] + dtp*u[1] + u[5] - u[3] - u[4])))/((1 + ^(u[2],2))*(3*dtdtp*dtp*^(u[1],2)*tauB*tauS - 3*^(dtp,2)*u[1]*tauB*tauS*^(u[2],2) + 3*dtdtp*dtp*^(u[1],2)*tauB*tauS*^(u[2],2) + 3*dmp*u[6]*tauB*tauS*(dtdmp*u[6] + dtdtp*u[1] + (-dtp + dtdmp*u[6] + dtdtp*u[1])*^(u[2],2)) - 4*dtdtp*u[1]*tauB*^(u[2],2)*visc - 3*dtdtp*u[1]*tauS*^(u[2],2)*zeta + 3*dtdtp*u[1]*tauB*tauS*u[5] - 3*dtp*tauB*tauS*^(u[2],2)*u[5] + 3*dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[5] - 3*dtdtp*u[1]*tauB*tauS*u[3] + 3*dtp*tauB*tauS*^(u[2],2)*u[3] + dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[3] + tauB*tauS*(3*dtp*^(u[2],2) + dtdtp*u[1]*(-3 + ^(u[2],2)))*u[4] + dtdmp*u[6]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))
            
            ,(-3*dtp*u[7]*(dtdmp*u[6] + dtdtp*u[1])*tauB*tauDiff*tauS*u[2])/(sqrt(1 + ^(u[2],2))*(3*dtdtp*dtp*^(u[1],2)*tauB*tauS - 3*^(dtp,2)*u[1]*tauB*tauS*^(u[2],2) + 3*dtdtp*dtp*^(u[1],2)*tauB*tauS*^(u[2],2) + 3*dmp*u[6]*tauB*tauS*(dtdmp*u[6] + dtdtp*u[1] + (-dtp + dtdmp*u[6] + dtdtp*u[1])*^(u[2],2)) - 4*dtdtp*u[1]*tauB*^(u[2],2)*visc - 3*dtdtp*u[1]*tauS*^(u[2],2)*zeta + 3*dtdtp*u[1]*tauB*tauS*u[5] - 3*dtp*tauB*tauS*^(u[2],2)*u[5] + 3*dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[5] - 3*dtdtp*u[1]*tauB*tauS*u[3] + 3*dtp*tauB*tauS*^(u[2],2)*u[3] + dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[3] + tauB*tauS*(3*dtp*^(u[2],2) + dtdtp*u[1]*(-3 + ^(u[2],2)))*u[4] + dtdmp*u[6]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))
            
            ,((1 + 2*^(u[2],2))*(dmp*u[6] + dtp*u[1] + u[5] - u[3] - u[4]))/sqrt(1 + ^(u[2],2))
            
            ,2*u[2]*(dmp*u[6] + dtp*u[1] + u[5] - u[3] - u[4])
            
            ,(-2*(visc - 2*tauS*u[3]))/(3*^(R,2))
            
            ,(-2*(visc - 2*tauS*u[4]))/(3*^(tau,2))
            
            ,zeta
            
            ,(-3*dtn*^(dtp,2)*^(u[1],2)*tauB*tauS + 3*dtdtp*dtp*n*^(u[1],2)*tauB*tauS + 3*^(dtp,2)*u[7]*u[1]*tauB*tauS*u[2] - 3*dtdtp*dtp*u[7]*^(u[1],2)*tauB*tauS*u[2] - 3*dtn*^(dtp,2)*^(u[1],2)*tauB*tauS*^(u[2],2) + 3*dtdtp*dtp*n*^(u[1],2)*tauB*tauS*^(u[2],2) - 3*^(dmp,2)*dtn*^(u[6],2)*tauB*tauS*(1 + ^(u[2],2)) + 4*dtdtp*u[7]*u[1]*tauB*u[2]*visc + 3*dtdtp*u[7]*u[1]*tauS*u[2]*zeta - 6*dtn*dtp*u[1]*tauB*tauS*u[5] + 3*dtdtp*n*u[1]*tauB*tauS*u[5] + 3*dtp*u[7]*tauB*tauS*u[2]*u[5] - 3*dtdtp*u[7]*u[1]*tauB*tauS*u[2]*u[5] - 6*dtn*dtp*u[1]*tauB*tauS*^(u[2],2)*u[5] + 3*dtdtp*n*u[1]*tauB*tauS*^(u[2],2)*u[5] - 3*dtn*tauB*tauS*^(u[5],2) - 3*dtn*tauB*tauS*^(u[2],2)*^(u[5],2) + 6*dtn*dtp*u[1]*tauB*tauS*u[3] - 3*dtdtp*n*u[1]*tauB*tauS*u[3] - 3*dtp*u[7]*tauB*tauS*u[2]*u[3] - dtdtp*u[7]*u[1]*tauB*tauS*u[2]*u[3] + 6*dtn*dtp*u[1]*tauB*tauS*^(u[2],2)*u[3] - 3*dtdtp*n*u[1]*tauB*tauS*^(u[2],2)*u[3] + 6*dtn*tauB*tauS*u[5]*u[3] + 6*dtn*tauB*tauS*^(u[2],2)*u[5]*u[3] - 3*dtn*tauB*tauS*^(u[3],2) - 3*dtn*tauB*tauS*^(u[2],2)*^(u[3],2) + 3*dmp*u[6]*tauB*tauS*(dtp*u[7]*u[2] + dtdmp*u[6]*(n - u[7]*u[2] + n*^(u[2],2)) + dtdtp*u[1]*(n - u[7]*u[2] + n*^(u[2],2)) - 2*dtn*(1 + ^(u[2],2))*(dtp*u[1] + u[5] - u[3] - u[4])) + tauB*tauS*(-3*dtp*u[7]*u[2] - dtdtp*u[1]*(u[7]*u[2] + 3*n*(1 + ^(u[2],2))) + 6*dtn*(1 + ^(u[2],2))*(dtp*u[1] + u[5] - u[3]))*u[4] - 3*dtn*tauB*tauS*(1 + ^(u[2],2))*^(u[4],2) + dtdmp*u[6]*(3*dtp*u[1]*tauB*tauS*(n - u[7]*u[2] + n*^(u[2],2)) + 3*n*tauB*tauS*(1 + ^(u[2],2))*(u[5] - u[3] - u[4]) + u[7]*u[2]*(4*tauB*visc + 3*tauS*zeta - tauB*tauS*(3*u[5] + u[3] + u[4]))))/((1 + ^(u[2],2))*(3*dtdtp*dtp*^(u[1],2)*tauB*tauS - 3*^(dtp,2)*u[1]*tauB*tauS*^(u[2],2) + 3*dtdtp*dtp*^(u[1],2)*tauB*tauS*^(u[2],2) + 3*dmp*u[6]*tauB*tauS*(dtdmp*u[6] + dtdtp*u[1] + (-dtp + dtdmp*u[6] + dtdtp*u[1])*^(u[2],2)) - 4*dtdtp*u[1]*tauB*^(u[2],2)*visc - 3*dtdtp*u[1]*tauS*^(u[2],2)*zeta + 3*dtdtp*u[1]*tauB*tauS*u[5] - 3*dtp*tauB*tauS*^(u[2],2)*u[5] + 3*dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[5] - 3*dtdtp*u[1]*tauB*tauS*u[3] + 3*dtp*tauB*tauS*^(u[2],2)*u[3] + dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[3] + tauB*tauS*(3*dtp*^(u[2],2) + dtdtp*u[1]*(-3 + ^(u[2],2)))*u[4] + dtdmp*u[6]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))
            
            ,-((u[7]*tauDiff*^(u[2],2)*(-3*^(dtp,2)*u[1]*tauB*tauS + 3*dmp*u[6]*(-dtp + dtdmp*u[6] + dtdtp*u[1])*tauB*tauS + 3*dtp*tauB*tauS*(dtdmp*u[6]*u[1] + dtdtp*^(u[1],2) - u[5] + u[3] + u[4]) - (dtdmp*u[6] + dtdtp*u[1])*(4*tauB*visc + 3*tauS*zeta - tauB*tauS*(3*u[5] + u[3] + u[4]))))/(sqrt(1 + ^(u[2],2))*(3*dtdtp*dtp*^(u[1],2)*tauB*tauS - 3*^(dtp,2)*u[1]*tauB*tauS*^(u[2],2) + 3*dtdtp*dtp*^(u[1],2)*tauB*tauS*^(u[2],2) + 3*dmp*u[6]*tauB*tauS*(dtdmp*u[6] + dtdtp*u[1] + (-dtp + dtdmp*u[6] + dtdtp*u[1])*^(u[2],2)) - 4*dtdtp*u[1]*tauB*^(u[2],2)*visc - 3*dtdtp*u[1]*tauS*^(u[2],2)*zeta + 3*dtdtp*u[1]*tauB*tauS*u[5] - 3*dtp*tauB*tauS*^(u[2],2)*u[5] + 3*dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[5] - 3*dtdtp*u[1]*tauB*tauS*u[3] + 3*dtp*tauB*tauS*^(u[2],2)*u[3] + dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[3] + tauB*tauS*(3*dtp*^(u[2],2) + dtdtp*u[1]*(-3 + ^(u[2],2)))*u[4] + dtdmp*u[6]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4])))))))
            
            ,-(u[2]*sqrt(1 + ^(u[2],2)))
            
            ,-1 - ^(u[2],2)
            
            ,(tauS*u[2])/^(R,2)
            
            ,0
            
            ,0
            
            ,(3*tauB*tauS*(dtdtp*u[1]*(u[7] + n*u[2] + n*^(u[2],3)) + dtdmp*u[6]*(u[7] + n*u[2]*(1 + ^(u[2],2))) + dtn*(u[2] + ^(u[2],3))*(-(dmp*u[6]) - dtp*u[1] - u[5] + u[3] + u[4])))/((1 + ^(u[2],2))*(3*dtdtp*dtp*^(u[1],2)*tauB*tauS - 3*^(dtp,2)*u[1]*tauB*tauS*^(u[2],2) + 3*dtdtp*dtp*^(u[1],2)*tauB*tauS*^(u[2],2) + 3*dmp*u[6]*tauB*tauS*(dtdmp*u[6] + dtdtp*u[1] + (-dtp + dtdmp*u[6] + dtdtp*u[1])*^(u[2],2)) - 4*dtdtp*u[1]*tauB*^(u[2],2)*visc - 3*dtdtp*u[1]*tauS*^(u[2],2)*zeta + 3*dtdtp*u[1]*tauB*tauS*u[5] - 3*dtp*tauB*tauS*^(u[2],2)*u[5] + 3*dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[5] - 3*dtdtp*u[1]*tauB*tauS*u[3] + 3*dtp*tauB*tauS*^(u[2],2)*u[3] + dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[3] + tauB*tauS*(3*dtp*^(u[2],2) + dtdtp*u[1]*(-3 + ^(u[2],2)))*u[4] + dtdmp*u[6]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))
            
            ,(3*u[7]*(dtdmp*u[6] + dtdtp*u[1])*tauB*tauDiff*tauS*u[2])/(sqrt(1 + ^(u[2],2))*(3*dtdtp*dtp*^(u[1],2)*tauB*tauS - 3*^(dtp,2)*u[1]*tauB*tauS*^(u[2],2) + 3*dtdtp*dtp*^(u[1],2)*tauB*tauS*^(u[2],2) + 3*dmp*u[6]*tauB*tauS*(dtdmp*u[6] + dtdtp*u[1] + (-dtp + dtdmp*u[6] + dtdtp*u[1])*^(u[2],2)) - 4*dtdtp*u[1]*tauB*^(u[2],2)*visc - 3*dtdtp*u[1]*tauS*^(u[2],2)*zeta + 3*dtdtp*u[1]*tauB*tauS*u[5] - 3*dtp*tauB*tauS*^(u[2],2)*u[5] + 3*dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[5] - 3*dtdtp*u[1]*tauB*tauS*u[3] + 3*dtp*tauB*tauS*^(u[2],2)*u[3] + dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[3] + tauB*tauS*(3*dtp*^(u[2],2) + dtdtp*u[1]*(-3 + ^(u[2],2)))*u[4] + dtdmp*u[6]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))
            
            ,-(u[2]*sqrt(1 + ^(u[2],2)))
            
            ,-1 - ^(u[2],2)
            
            ,0
            
            ,(tauS*u[2])/^(tau,2)
            
            ,0
            
            ,(3*tauB*tauS*(dtdtp*u[1]*(u[7] + n*u[2] + n*^(u[2],3)) + dtdmp*u[6]*(u[7] + n*u[2]*(1 + ^(u[2],2))) + dtn*(u[2] + ^(u[2],3))*(-(dmp*u[6]) - dtp*u[1] - u[5] + u[3] + u[4])))/((1 + ^(u[2],2))*(3*dtdtp*dtp*^(u[1],2)*tauB*tauS - 3*^(dtp,2)*u[1]*tauB*tauS*^(u[2],2) + 3*dtdtp*dtp*^(u[1],2)*tauB*tauS*^(u[2],2) + 3*dmp*u[6]*tauB*tauS*(dtdmp*u[6] + dtdtp*u[1] + (-dtp + dtdmp*u[6] + dtdtp*u[1])*^(u[2],2)) - 4*dtdtp*u[1]*tauB*^(u[2],2)*visc - 3*dtdtp*u[1]*tauS*^(u[2],2)*zeta + 3*dtdtp*u[1]*tauB*tauS*u[5] - 3*dtp*tauB*tauS*^(u[2],2)*u[5] + 3*dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[5] - 3*dtdtp*u[1]*tauB*tauS*u[3] + 3*dtp*tauB*tauS*^(u[2],2)*u[3] + dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[3] + tauB*tauS*(3*dtp*^(u[2],2) + dtdtp*u[1]*(-3 + ^(u[2],2)))*u[4] + dtdmp*u[6]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))
            
            ,(3*u[7]*(dtdmp*u[6] + dtdtp*u[1])*tauB*tauDiff*tauS*u[2])/(sqrt(1 + ^(u[2],2))*(3*dtdtp*dtp*^(u[1],2)*tauB*tauS - 3*^(dtp,2)*u[1]*tauB*tauS*^(u[2],2) + 3*dtdtp*dtp*^(u[1],2)*tauB*tauS*^(u[2],2) + 3*dmp*u[6]*tauB*tauS*(dtdmp*u[6] + dtdtp*u[1] + (-dtp + dtdmp*u[6] + dtdtp*u[1])*^(u[2],2)) - 4*dtdtp*u[1]*tauB*^(u[2],2)*visc - 3*dtdtp*u[1]*tauS*^(u[2],2)*zeta + 3*dtdtp*u[1]*tauB*tauS*u[5] - 3*dtp*tauB*tauS*^(u[2],2)*u[5] + 3*dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[5] - 3*dtdtp*u[1]*tauB*tauS*u[3] + 3*dtp*tauB*tauS*^(u[2],2)*u[3] + dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[3] + tauB*tauS*(3*dtp*^(u[2],2) + dtdtp*u[1]*(-3 + ^(u[2],2)))*u[4] + dtdmp*u[6]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))
            
            ,u[2]*sqrt(1 + ^(u[2],2))
            
            ,1 + ^(u[2],2)
            
            ,0
            
            ,0
            
            ,tauB*u[2]
            
            ,(-3*tauB*tauS*(dtdtp*u[1]*(u[7] + n*u[2] + n*^(u[2],3)) + dtdmp*u[6]*(u[7] + n*u[2]*(1 + ^(u[2],2))) + dtn*(u[2] + ^(u[2],3))*(-(dmp*u[6]) - dtp*u[1] - u[5] + u[3] + u[4])))/((1 + ^(u[2],2))*(3*dtdtp*dtp*^(u[1],2)*tauB*tauS - 3*^(dtp,2)*u[1]*tauB*tauS*^(u[2],2) + 3*dtdtp*dtp*^(u[1],2)*tauB*tauS*^(u[2],2) + 3*dmp*u[6]*tauB*tauS*(dtdmp*u[6] + dtdtp*u[1] + (-dtp + dtdmp*u[6] + dtdtp*u[1])*^(u[2],2)) - 4*dtdtp*u[1]*tauB*^(u[2],2)*visc - 3*dtdtp*u[1]*tauS*^(u[2],2)*zeta + 3*dtdtp*u[1]*tauB*tauS*u[5] - 3*dtp*tauB*tauS*^(u[2],2)*u[5] + 3*dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[5] - 3*dtdtp*u[1]*tauB*tauS*u[3] + 3*dtp*tauB*tauS*^(u[2],2)*u[3] + dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[3] + tauB*tauS*(3*dtp*^(u[2],2) + dtdtp*u[1]*(-3 + ^(u[2],2)))*u[4] + dtdmp*u[6]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))
            
            ,(-3*u[7]*(dtdmp*u[6] + dtdtp*u[1])*tauB*tauDiff*tauS*u[2])/(sqrt(1 + ^(u[2],2))*(3*dtdtp*dtp*^(u[1],2)*tauB*tauS - 3*^(dtp,2)*u[1]*tauB*tauS*^(u[2],2) + 3*dtdtp*dtp*^(u[1],2)*tauB*tauS*^(u[2],2) + 3*dmp*u[6]*tauB*tauS*(dtdmp*u[6] + dtdtp*u[1] + (-dtp + dtdmp*u[6] + dtdtp*u[1])*^(u[2],2)) - 4*dtdtp*u[1]*tauB*^(u[2],2)*visc - 3*dtdtp*u[1]*tauS*^(u[2],2)*zeta + 3*dtdtp*u[1]*tauB*tauS*u[5] - 3*dtp*tauB*tauS*^(u[2],2)*u[5] + 3*dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[5] - 3*dtdtp*u[1]*tauB*tauS*u[3] + 3*dtp*tauB*tauS*^(u[2],2)*u[3] + dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[3] + tauB*tauS*(3*dtp*^(u[2],2) + dtdtp*u[1]*(-3 + ^(u[2],2)))*u[4] + dtdmp*u[6]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))
            
            ,(dmp + dmdmp*u[6] + dtdmp*u[1])*u[2]*sqrt(1 + ^(u[2],2))
            
            ,dmp + (dmp + dmdmp*u[6] + dtdmp*u[1])*^(u[2],2)
            
            ,0
            
            ,0
            
            ,0
            
            ,(3*^(dmp,2)*dtn*u[6]*tauB*tauS*u[2]*^(1 + ^(u[2],2),2) + 3*dmp*tauB*tauS*(1 + ^(u[2],2))*(-(dtdmp*u[6]*u[7]) - dtdtp*u[7]*u[1] - dmn*dtp*u[6]*^(u[2],3) - dmdmp*dtn*^(u[6],2)*u[2]*(1 + ^(u[2],2)) + dtdtp*(dmn*u[6] - n)*u[1]*u[2]*(1 + ^(u[2],2)) + dtdmp*u[6]*(dmn*u[6] - n - dtn*u[1])*u[2]*(1 + ^(u[2],2)) + dtn*u[2]*(1 + ^(u[2],2))*(dtp*u[1] + u[5] - u[3] - u[4])) + u[2]*(dmn*(1 + ^(u[2],2))*(3*dtp*tauB*tauS*^(u[2],2)*(-(dtp*u[1]) - u[5] + u[3] + u[4]) + dtdmp*u[6]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4])))) + dtdtp*u[1]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))) - (dmdmp*u[6] + dtdmp*u[1])*(-3*dtp*tauB*tauS*u[2]*(u[7] + n*u[2]*(1 + ^(u[2],2))) + dtn*(1 + ^(u[2],2))*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4])))))))/((1 + ^(u[2],2))*(3*dtdtp*dtp*^(u[1],2)*tauB*tauS - 3*^(dtp,2)*u[1]*tauB*tauS*^(u[2],2) + 3*dtdtp*dtp*^(u[1],2)*tauB*tauS*^(u[2],2) + 3*dmp*u[6]*tauB*tauS*(dtdmp*u[6] + dtdtp*u[1] + (-dtp + dtdmp*u[6] + dtdtp*u[1])*^(u[2],2)) - 4*dtdtp*u[1]*tauB*^(u[2],2)*visc - 3*dtdtp*u[1]*tauS*^(u[2],2)*zeta + 3*dtdtp*u[1]*tauB*tauS*u[5] - 3*dtp*tauB*tauS*^(u[2],2)*u[5] + 3*dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[5] - 3*dtdtp*u[1]*tauB*tauS*u[3] + 3*dtp*tauB*tauS*^(u[2],2)*u[3] + dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[3] + tauB*tauS*(3*dtp*^(u[2],2) + dtdtp*u[1]*(-3 + ^(u[2],2)))*u[4] + dtdmp*u[6]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))
            
            ,(-3*dmp*dtdmp*u[6]*u[7]*tauB*tauDiff*tauS*u[2] - 3*dmp*dtdtp*u[7]*u[1]*tauB*tauDiff*tauS*u[2] - 3*dmp*dtdmp*u[6]*u[7]*tauB*tauDiff*tauS*^(u[2],3) + 3*dmdmp*dtp*u[6]*u[7]*tauB*tauDiff*tauS*^(u[2],3) - 3*dmp*dtdtp*u[7]*u[1]*tauB*tauDiff*tauS*^(u[2],3) + 3*dtdmp*dtp*u[7]*u[1]*tauB*tauDiff*tauS*^(u[2],3) + Ds*n*sqrt(1 + ^(u[2],2))*(3*(dtdmp*u[6] + dtdtp*u[1])*tauB*tauS*(dmp*u[6] + dtp*u[1] + u[5] - u[3] - u[4]) + ^(u[2],2)*(-3*^(dtp,2)*u[1]*tauB*tauS + 3*dmp*u[6]*(-dtp + dtdmp*u[6] + dtdtp*u[1])*tauB*tauS + 3*dtp*tauB*tauS*(dtdmp*u[6]*u[1] + dtdtp*^(u[1],2) - u[5] + u[3] + u[4]) - (dtdmp*u[6] + dtdtp*u[1])*(4*tauB*visc + 3*tauS*zeta - tauB*tauS*(3*u[5] + u[3] + u[4])))) + Ds*n*^(u[2],2)*sqrt(1 + ^(u[2],2))*(3*(dtdmp*u[6] + dtdtp*u[1])*tauB*tauS*(dmp*u[6] + dtp*u[1] + u[5] - u[3] - u[4]) + ^(u[2],2)*(-3*^(dtp,2)*u[1]*tauB*tauS + 3*dmp*u[6]*(-dtp + dtdmp*u[6] + dtdtp*u[1])*tauB*tauS + 3*dtp*tauB*tauS*(dtdmp*u[6]*u[1] + dtdtp*^(u[1],2) - u[5] + u[3] + u[4]) - (dtdmp*u[6] + dtdtp*u[1])*(4*tauB*visc + 3*tauS*zeta - tauB*tauS*(3*u[5] + u[3] + u[4])))))/(sqrt(1 + ^(u[2],2))*(3*(dtdmp*u[6] + dtdtp*u[1])*tauB*tauS*(dmp*u[6] + dtp*u[1] + u[5] - u[3] - u[4]) + ^(u[2],2)*(-3*^(dtp,2)*u[1]*tauB*tauS + 3*dmp*u[6]*(-dtp + dtdmp*u[6] + dtdtp*u[1])*tauB*tauS + 3*dtp*tauB*tauS*(dtdmp*u[6]*u[1] + dtdtp*^(u[1],2) - u[5] + u[3] + u[4]) - (dtdmp*u[6] + dtdtp*u[1])*(4*tauB*visc + 3*tauS*zeta - tauB*tauS*(3*u[5] + u[3] + u[4])))))
            
            ,0
            
            ,0
            
            ,0
            
            ,0
            
            ,0
            
            ,1
            
            ,0
        )
            #########################################################################################################################################################################
        
        source=SVector{7}(
            (dmp*u[6]*tau*u[2]*sqrt(1 + ^(u[2],2)) + dtp*u[1]*tau*u[2]*sqrt(1 + ^(u[2],2)) + tau*u[2]*sqrt(1 + ^(u[2],2))*u[5] - tau*u[2]*sqrt(1 + ^(u[2],2))*u[3] + R*^(u[2],2)*(dmp*u[6] + dtp*u[1] + u[5] - u[3] - u[4]) - tau*u[2]*sqrt(1 + ^(u[2],2))*u[4] + R*(dmp*u[6] + dtp*u[1] + u[5] + u[4]))/(R*tau)
    
            ,(dmp*u[6]*u[2]*(tau*u[2] + R*sqrt(1 + ^(u[2],2))) + dtp*u[1]*u[2]*(tau*u[2] + R*sqrt(1 + ^(u[2],2))) + tau*^(u[2],2)*u[5] + R*u[2]*sqrt(1 + ^(u[2],2))*u[5] - R*u[2]*sqrt(1 + ^(u[2],2))*u[3] - R*u[2]*sqrt(1 + ^(u[2],2))*u[4] - tau*((2 + ^(u[2],2))*u[3] + u[4] + ^(u[2],2)*u[4]))/(R*tau)
    
            ,(3*R*u[3] - (2*R*sqrt(1 + ^(u[2],2))*(visc - 2*tauS*u[3]))/tau + 4*u[2]*(visc + tauS*u[3]))/(3*^(R,3))
    
            ,(3*tau*u[4] - (2*tau*u[2]*(visc - 2*tauS*u[4]))/R + 4*sqrt(1 + ^(u[2],2))*(visc + tauS*u[4]))/(3*^(tau,3))
    
            ,(u[2]*zeta)/R + (sqrt(1 + ^(u[2],2))*zeta)/tau + u[5]
    
            ,-((-3*dtdtp*dtp*n*R*^(u[1],2)*tauB*tauS - 3*^(dtp,2)*u[7]*R*u[1]*tauB*tauS*u[2] - 3*dtdtp*dtp*u[7]*R*^(u[1],2)*tauB*tauS*u[2] - 9*dtdtp*dtp*n*R*^(u[1],2)*tauB*tauS*^(u[2],2) - 6*dtdtp*dtp*u[7]*R*^(u[1],2)*tauB*tauS*^(u[2],3) - 9*dtdtp*dtp*n*R*^(u[1],2)*tauB*tauS*^(u[2],4) + 3*^(dtp,2)*u[7]*R*u[1]*tauB*tauS*^(u[2],5) - 3*dtdtp*dtp*u[7]*R*^(u[1],2)*tauB*tauS*^(u[2],5) - 3*dtdtp*dtp*n*R*^(u[1],2)*tauB*tauS*^(u[2],6) - 3*dtdtp*dtp*u[7]*^(u[1],2)*tau*tauB*tauS*sqrt(1 + ^(u[2],2)) + 3*dtn*^(dtp,2)*^(u[1],2)*tau*tauB*tauS*u[2]*sqrt(1 + ^(u[2],2)) - 3*dtdtp*dtp*n*^(u[1],2)*tau*tauB*tauS*u[2]*sqrt(1 + ^(u[2],2)) - 6*dtdtp*dtp*u[7]*^(u[1],2)*tau*tauB*tauS*^(u[2],2)*sqrt(1 + ^(u[2],2)) + 6*dtn*^(dtp,2)*^(u[1],2)*tau*tauB*tauS*^(u[2],3)*sqrt(1 + ^(u[2],2)) - 6*dtdtp*dtp*n*^(u[1],2)*tau*tauB*tauS*^(u[2],3)*sqrt(1 + ^(u[2],2)) + 3*^(dtp,2)*u[7]*u[1]*tau*tauB*tauS*^(u[2],4)*sqrt(1 + ^(u[2],2)) - 3*dtdtp*dtp*u[7]*^(u[1],2)*tau*tauB*tauS*^(u[2],4)*sqrt(1 + ^(u[2],2)) + 3*dtn*^(dtp,2)*^(u[1],2)*tau*tauB*tauS*^(u[2],5)*sqrt(1 + ^(u[2],2)) - 3*dtdtp*dtp*n*^(u[1],2)*tau*tauB*tauS*^(u[2],5)*sqrt(1 + ^(u[2],2)) + 3*^(dmp,2)*dtn*^(u[6],2)*tauB*tauS*^(1 + ^(u[2],2),2)*(R + R*^(u[2],2) + tau*u[2]*sqrt(1 + ^(u[2],2))) + 2*dtdtp*u[7]*R*u[1]*tauB*u[2]*visc + 6*dtdtp*n*R*u[1]*tauB*^(u[2],2)*visc + 6*dtdtp*u[7]*R*u[1]*tauB*^(u[2],3)*visc + 12*dtdtp*n*R*u[1]*tauB*^(u[2],4)*visc + 4*dtdtp*u[7]*R*u[1]*tauB*^(u[2],5)*visc + 6*dtdtp*n*R*u[1]*tauB*^(u[2],6)*visc + 6*dtdtp*u[7]*u[1]*tau*tauB*^(u[2],2)*sqrt(1 + ^(u[2],2))*visc - 6*dtn*dtp*u[1]*tau*tauB*^(u[2],3)*sqrt(1 + ^(u[2],2))*visc + 6*dtdtp*n*u[1]*tau*tauB*^(u[2],3)*sqrt(1 + ^(u[2],2))*visc + 4*dtdtp*u[7]*u[1]*tau*tauB*^(u[2],4)*sqrt(1 + ^(u[2],2))*visc - 6*dtn*dtp*u[1]*tau*tauB*^(u[2],5)*sqrt(1 + ^(u[2],2))*visc + 6*dtdtp*n*u[1]*tau*tauB*^(u[2],5)*sqrt(1 + ^(u[2],2))*visc - 3*dtdtp*u[7]*R*u[1]*tauS*u[2]*zeta + 3*dtdtp*u[7]*R*u[1]*tauS*^(u[2],5)*zeta + 3*dtdtp*u[7]*u[1]*tau*tauS*^(u[2],4)*sqrt(1 + ^(u[2],2))*zeta - 3*dtdtp*n*R*u[1]*tauB*tauS*u[5] - 3*dtp*u[7]*R*tauB*tauS*u[2]*u[5] - 3*dtdtp*u[7]*R*u[1]*tauB*tauS*u[2]*u[5] - 9*dtdtp*n*R*u[1]*tauB*tauS*^(u[2],2)*u[5] - 6*dtdtp*u[7]*R*u[1]*tauB*tauS*^(u[2],3)*u[5] - 9*dtdtp*n*R*u[1]*tauB*tauS*^(u[2],4)*u[5] + 3*dtp*u[7]*R*tauB*tauS*^(u[2],5)*u[5] - 3*dtdtp*u[7]*R*u[1]*tauB*tauS*^(u[2],5)*u[5] - 3*dtdtp*n*R*u[1]*tauB*tauS*^(u[2],6)*u[5] - 3*dtdtp*u[7]*u[1]*tau*tauB*tauS*sqrt(1 + ^(u[2],2))*u[5] - 3*dtdtp*u[7]*R*u[1]*tau*tauS*u[2]*sqrt(1 + ^(u[2],2))*u[5] + 6*dtn*dtp*u[1]*tau*tauB*tauS*u[2]*sqrt(1 + ^(u[2],2))*u[5] - 3*dtdtp*n*u[1]*tau*tauB*tauS*u[2]*sqrt(1 + ^(u[2],2))*u[5] + 3*dtn*dtp*R*u[1]*tau*tauS*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[5] - 3*dtdtp*n*R*u[1]*tau*tauS*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[5] - 6*dtdtp*u[7]*u[1]*tau*tauB*tauS*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[5] + 12*dtn*dtp*u[1]*tau*tauB*tauS*^(u[2],3)*sqrt(1 + ^(u[2],2))*u[5] - 6*dtdtp*n*u[1]*tau*tauB*tauS*^(u[2],3)*sqrt(1 + ^(u[2],2))*u[5] + 3*dtn*dtp*R*u[1]*tau*tauS*^(u[2],4)*sqrt(1 + ^(u[2],2))*u[5] - 3*dtdtp*n*R*u[1]*tau*tauS*^(u[2],4)*sqrt(1 + ^(u[2],2))*u[5] + 3*dtp*u[7]*tau*tauB*tauS*^(u[2],4)*sqrt(1 + ^(u[2],2))*u[5] - 3*dtdtp*u[7]*u[1]*tau*tauB*tauS*^(u[2],4)*sqrt(1 + ^(u[2],2))*u[5] + 6*dtn*dtp*u[1]*tau*tauB*tauS*^(u[2],5)*sqrt(1 + ^(u[2],2))*u[5] - 3*dtdtp*n*u[1]*tau*tauB*tauS*^(u[2],5)*sqrt(1 + ^(u[2],2))*u[5] - 6*dtn*tau*tauB*^(u[2],3)*sqrt(1 + ^(u[2],2))*visc*u[5] - 6*dtn*tau*tauB*^(u[2],5)*sqrt(1 + ^(u[2],2))*visc*u[5] + 3*dtn*tau*tauB*tauS*u[2]*sqrt(1 + ^(u[2],2))*^(u[5],2) + 3*dtn*R*tau*tauS*^(u[2],2)*sqrt(1 + ^(u[2],2))*^(u[5],2) + 6*dtn*tau*tauB*tauS*^(u[2],3)*sqrt(1 + ^(u[2],2))*^(u[5],2) + 3*dtn*R*tau*tauS*^(u[2],4)*sqrt(1 + ^(u[2],2))*^(u[5],2) + 3*dtn*tau*tauB*tauS*^(u[2],5)*sqrt(1 + ^(u[2],2))*^(u[5],2) + 3*dtdtp*n*R*u[1]*tauB*tauS*u[3] + 4*dtdtp*u[7]*R*u[1]*tauB*tauS*u[2]*u[3] - 3*dtp*n*R*tauB*tauS*^(u[2],2)*u[3] + 6*dtdtp*n*R*u[1]*tauB*tauS*^(u[2],2)*u[3] - 3*dtp*u[7]*R*tauB*tauS*^(u[2],3)*u[3] + 3*dtdtp*u[7]*R*u[1]*tauB*tauS*^(u[2],3)*u[3] - 6*dtp*n*R*tauB*tauS*^(u[2],4)*u[3] + 3*dtdtp*n*R*u[1]*tauB*tauS*^(u[2],4)*u[3] - 3*dtp*u[7]*R*tauB*tauS*^(u[2],5)*u[3] - dtdtp*u[7]*R*u[1]*tauB*tauS*^(u[2],5)*u[3] - 3*dtp*n*R*tauB*tauS*^(u[2],6)*u[3] - 3*dtdtp*u[7]*u[1]*tau*tauB*tauS*sqrt(1 + ^(u[2],2))*u[3] + 3*dtdtp*u[7]*R*u[1]*tau*tauB*u[2]*sqrt(1 + ^(u[2],2))*u[3] + 6*dtn*dtp*u[1]*tau*tauB*tauS*u[2]*sqrt(1 + ^(u[2],2))*u[3] - 3*dtdtp*n*u[1]*tau*tauB*tauS*u[2]*sqrt(1 + ^(u[2],2))*u[3] - 3*dtn*dtp*R*u[1]*tau*tauB*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[3] + 3*dtdtp*n*R*u[1]*tau*tauB*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[3] - 6*dtp*u[7]*tau*tauB*tauS*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[3] - 6*dtp*n*tau*tauB*tauS*^(u[2],3)*sqrt(1 + ^(u[2],2))*u[3] + 12*dtn*dtp*u[1]*tau*tauB*tauS*^(u[2],3)*sqrt(1 + ^(u[2],2))*u[3] - 6*dtdtp*n*u[1]*tau*tauB*tauS*^(u[2],3)*sqrt(1 + ^(u[2],2))*u[3] - 3*dtn*dtp*R*u[1]*tau*tauB*^(u[2],4)*sqrt(1 + ^(u[2],2))*u[3] + 3*dtdtp*n*R*u[1]*tau*tauB*^(u[2],4)*sqrt(1 + ^(u[2],2))*u[3] - 3*dtp*u[7]*tau*tauB*tauS*^(u[2],4)*sqrt(1 + ^(u[2],2))*u[3] - dtdtp*u[7]*u[1]*tau*tauB*tauS*^(u[2],4)*sqrt(1 + ^(u[2],2))*u[3] - 6*dtp*n*tau*tauB*tauS*^(u[2],5)*sqrt(1 + ^(u[2],2))*u[3] + 6*dtn*dtp*u[1]*tau*tauB*tauS*^(u[2],5)*sqrt(1 + ^(u[2],2))*u[3] - 3*dtdtp*n*u[1]*tau*tauB*tauS*^(u[2],5)*sqrt(1 + ^(u[2],2))*u[3] - 2*dtn*tau*tauB*^(u[2],3)*sqrt(1 + ^(u[2],2))*visc*u[3] - 2*dtn*tau*tauB*^(u[2],5)*sqrt(1 + ^(u[2],2))*visc*u[3] - 6*dtn*tau*tauS*^(u[2],3)*sqrt(1 + ^(u[2],2))*zeta*u[3] - 6*dtn*tau*tauS*^(u[2],5)*sqrt(1 + ^(u[2],2))*zeta*u[3] + 6*dtn*tau*tauB*tauS*u[2]*sqrt(1 + ^(u[2],2))*u[5]*u[3] - 3*dtn*R*tau*tauB*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[5]*u[3] - 3*dtn*R*tau*tauS*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[5]*u[3] + 12*dtn*tau*tauB*tauS*^(u[2],3)*sqrt(1 + ^(u[2],2))*u[5]*u[3] - 3*dtn*R*tau*tauB*^(u[2],4)*sqrt(1 + ^(u[2],2))*u[5]*u[3] - 3*dtn*R*tau*tauS*^(u[2],4)*sqrt(1 + ^(u[2],2))*u[5]*u[3] + 6*dtn*tau*tauB*tauS*^(u[2],5)*sqrt(1 + ^(u[2],2))*u[5]*u[3] - 9*dtn*tau*tauB*tauS*u[2]*sqrt(1 + ^(u[2],2))*^(u[3],2) + 3*dtn*R*tau*tauB*^(u[2],2)*sqrt(1 + ^(u[2],2))*^(u[3],2) - 10*dtn*tau*tauB*tauS*^(u[2],3)*sqrt(1 + ^(u[2],2))*^(u[3],2) + 3*dtn*R*tau*tauB*^(u[2],4)*sqrt(1 + ^(u[2],2))*^(u[3],2) - dtn*tau*tauB*tauS*^(u[2],5)*sqrt(1 + ^(u[2],2))*^(u[3],2) + 3*dtdtp*n*R*u[1]*tauB*tauS*u[4] - 3*dtp*u[7]*R*tauB*tauS*u[2]*u[4] + dtdtp*u[7]*R*u[1]*tauB*tauS*u[2]*u[4] - 6*dtp*n*R*tauB*tauS*^(u[2],2)*u[4] + 3*dtdtp*n*R*u[1]*tauB*tauS*^(u[2],2)*u[4] - 6*dtp*u[7]*R*tauB*tauS*^(u[2],3)*u[4] - 12*dtp*n*R*tauB*tauS*^(u[2],4)*u[4] - 3*dtdtp*n*R*u[1]*tauB*tauS*^(u[2],4)*u[4] - 3*dtp*u[7]*R*tauB*tauS*^(u[2],5)*u[4] - dtdtp*u[7]*R*u[1]*tauB*tauS*^(u[2],5)*u[4] - 6*dtp*n*R*tauB*tauS*^(u[2],6)*u[4] - 3*dtdtp*n*R*u[1]*tauB*tauS*^(u[2],6)*u[4] + 3*dtdtp*u[7]*R*u[1]*tau*tauB*u[2]*sqrt(1 + ^(u[2],2))*u[4] - 3*dtn*dtp*R*u[1]*tau*tauB*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[4] + 3*dtdtp*n*R*u[1]*tau*tauB*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[4] - 3*dtp*u[7]*tau*tauB*tauS*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[4] + 3*dtdtp*u[7]*u[1]*tau*tauB*tauS*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[4] - 3*dtp*n*tau*tauB*tauS*^(u[2],3)*sqrt(1 + ^(u[2],2))*u[4] - 3*dtn*dtp*R*u[1]*tau*tauB*^(u[2],4)*sqrt(1 + ^(u[2],2))*u[4] + 3*dtdtp*n*R*u[1]*tau*tauB*^(u[2],4)*sqrt(1 + ^(u[2],2))*u[4] - 3*dtp*u[7]*tau*tauB*tauS*^(u[2],4)*sqrt(1 + ^(u[2],2))*u[4] - dtdtp*u[7]*u[1]*tau*tauB*tauS*^(u[2],4)*sqrt(1 + ^(u[2],2))*u[4] - 3*dtp*n*tau*tauB*tauS*^(u[2],5)*sqrt(1 + ^(u[2],2))*u[4] + 2*dtn*tau*tauB*^(u[2],3)*sqrt(1 + ^(u[2],2))*visc*u[4] + 2*dtn*tau*tauB*^(u[2],5)*sqrt(1 + ^(u[2],2))*visc*u[4] - 3*dtn*tau*tauS*^(u[2],3)*sqrt(1 + ^(u[2],2))*zeta*u[4] - 3*dtn*tau*tauS*^(u[2],5)*sqrt(1 + ^(u[2],2))*zeta*u[4] - 3*dtn*R*tau*tauB*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[5]*u[4] - 3*dtn*R*tau*tauS*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[5]*u[4] - 3*dtn*R*tau*tauB*^(u[2],4)*sqrt(1 + ^(u[2],2))*u[5]*u[4] - 3*dtn*R*tau*tauS*^(u[2],4)*sqrt(1 + ^(u[2],2))*u[5]*u[4] - 12*dtn*tau*tauB*tauS*u[2]*sqrt(1 + ^(u[2],2))*u[3]*u[4] + 6*dtn*R*tau*tauB*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[3]*u[4] - 12*dtn*tau*tauB*tauS*^(u[2],3)*sqrt(1 + ^(u[2],2))*u[3]*u[4] + 6*dtn*R*tau*tauB*^(u[2],4)*sqrt(1 + ^(u[2],2))*u[3]*u[4] - 3*dtn*tau*tauB*tauS*u[2]*sqrt(1 + ^(u[2],2))*^(u[4],2) + 3*dtn*R*tau*tauB*^(u[2],2)*sqrt(1 + ^(u[2],2))*^(u[4],2) - 2*dtn*tau*tauB*tauS*^(u[2],3)*sqrt(1 + ^(u[2],2))*^(u[4],2) + 3*dtn*R*tau*tauB*^(u[2],4)*sqrt(1 + ^(u[2],2))*^(u[4],2) + dtn*tau*tauB*tauS*^(u[2],5)*sqrt(1 + ^(u[2],2))*^(u[4],2) + dtn*R*^(1 + ^(u[2],2),2)*(3*^(dtp,2)*^(u[1],2)*tauB*tauS*(1 + ^(u[2],2)) + 3*tauB*tauS*(u[5] - u[3] - u[4])*(u[5] + u[4]) - 3*tauS*^(u[2],2)*zeta*(u[3] + 2*u[4]) + tauB*^(u[2],2)*(-2*visc*(3*u[5] - u[3] + u[4]) + tauS*(3*^(u[5],2) + ^(u[3],2) + 6*u[5]*u[4] - ^(u[4],2))) + 3*dtp*u[1]*tauB*(2*tauS*u[5] - tauS*u[3] + 2*^(u[2],2)*(-visc + tauS*(u[5] + u[4])))) - 3*dmp*u[6]*(dtdmp*u[6]*tauB*tauS*^(1 + ^(u[2],2),2)*(u[7]*R*u[2] + u[7]*tau*sqrt(1 + ^(u[2],2)) + n*tau*u[2]*sqrt(1 + ^(u[2],2)) + n*R*(1 + ^(u[2],2))) + tauB*tauS*(dtp*u[7]*u[2]*(R - R*^(u[2],4) - tau*^(u[2],3)*sqrt(1 + ^(u[2],2))) + dtdtp*u[1]*^(1 + ^(u[2],2),2)*(u[7]*R*u[2] + u[7]*tau*sqrt(1 + ^(u[2],2)) + n*tau*u[2]*sqrt(1 + ^(u[2],2)) + n*R*(1 + ^(u[2],2)))) - dtn*(1 + ^(u[2],2))*(2*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2))*(R + R*^(u[2],2) + tau*u[2]*sqrt(1 + ^(u[2],2))) + 2*tau*tauB*u[2]*sqrt(1 + ^(u[2],2))*(tauS*(u[5] + u[3]) + ^(u[2],2)*(-visc + tauS*(u[5] + u[3]))) + R*(tau*tauS*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[5] + tauB*(2*tauS*u[5] - tauS*u[3] + 2*^(u[2],4)*(-visc + tauS*(u[5] + u[4])) - ^(u[2],2)*(2*visc + tauS*(-4*u[5] + u[3] - 2*u[4]) + tau*sqrt(1 + ^(u[2],2))*(u[3] + u[4])))))) - dtdmp*u[6]*(3*dtp*u[1]*tauB*tauS*^(1 + ^(u[2],2),2)*(u[7]*R*u[2] + u[7]*tau*sqrt(1 + ^(u[2],2)) + n*tau*u[2]*sqrt(1 + ^(u[2],2)) + n*R*(1 + ^(u[2],2))) - 3*n*(1 + ^(u[2],2))*(-(tau*tauB*u[2]*sqrt(1 + ^(u[2],2))*(tauS*(u[5] + u[3]) + ^(u[2],2)*(-2*visc + tauS*(u[5] + u[3])))) + R*(-(tau*tauS*^(u[2],2)*sqrt(1 + ^(u[2],2))*u[5]) + tauB*(tauS*(-u[5] + u[3] + u[4]) + ^(u[2],4)*(2*visc - tauS*(u[5] + u[4])) + ^(u[2],2)*(2*visc + tauS*(-2*u[5] + u[3]) + tau*sqrt(1 + ^(u[2],2))*(u[3] + u[4]))))) + u[7]*(tau*sqrt(1 + ^(u[2],2))*(-3*tauS*^(u[2],4)*zeta + tauB*(3*tauS*(u[5] + u[3]) - 3*^(u[2],2)*(2*visc - 2*tauS*u[5] + tauS*u[4]) + ^(u[2],4)*(-4*visc + tauS*(3*u[5] + u[3] + u[4])))) + R*(3*tauS*u[2]*(zeta - ^(u[2],4)*zeta + tau*sqrt(1 + ^(u[2],2))*u[5]) + tauB*(-3*^(u[2],3)*(2*visc + tauS*(-2*u[5] + u[3])) + ^(u[2],5)*(-4*visc + tauS*(3*u[5] + u[3] + u[4])) - u[2]*(2*visc + 3*tau*sqrt(1 + ^(u[2],2))*(u[3] + u[4]) + tauS*(-3*u[5] + 4*u[3] + u[4])))))))/(R*tau*^(1 + ^(u[2],2),1.5)*(3*dtdtp*dtp*^(u[1],2)*tauB*tauS - 3*^(dtp,2)*u[1]*tauB*tauS*^(u[2],2) + 3*dtdtp*dtp*^(u[1],2)*tauB*tauS*^(u[2],2) + 3*dmp*u[6]*tauB*tauS*(dtdmp*u[6] + dtdtp*u[1] + (-dtp + dtdmp*u[6] + dtdtp*u[1])*^(u[2],2)) - 4*dtdtp*u[1]*tauB*^(u[2],2)*visc - 3*dtdtp*u[1]*tauS*^(u[2],2)*zeta + 3*dtdtp*u[1]*tauB*tauS*u[5] - 3*dtp*tauB*tauS*^(u[2],2)*u[5] + 3*dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[5] - 3*dtdtp*u[1]*tauB*tauS*u[3] + 3*dtp*tauB*tauS*^(u[2],2)*u[3] + dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[3] + tauB*tauS*(3*dtp*^(u[2],2) + dtdtp*u[1]*(-3 + ^(u[2],2)))*u[4] + dtdmp*u[6]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4])))))))
    
            ,(u[7]*(3*^(dtp,2)*u[1]*tauB*tauS*^(u[2],2)*(tau*tauDiff*u[2] - R*tau*sqrt(1 + ^(u[2],2)) + R*tauDiff*sqrt(1 + ^(u[2],2))) + 3*dmp*u[6]*tauB*tauS*(R*(dtdmp*u[6] + dtdtp*u[1])*tau*^(1 + ^(u[2],2),1.5) + dtp*^(u[2],2)*(tau*tauDiff*u[2] - R*tau*sqrt(1 + ^(u[2],2)) + R*tauDiff*sqrt(1 + ^(u[2],2)))) + 3*dtp*tauB*tauS*(dtdmp*u[6]*R*u[1]*tau*^(1 + ^(u[2],2),1.5) + dtdtp*R*^(u[1],2)*tau*^(1 + ^(u[2],2),1.5) + ^(u[2],2)*(tau*tauDiff*u[2]*(u[5] + u[3]) + R*tauDiff*sqrt(1 + ^(u[2],2))*(u[5] + u[4]) + R*tau*sqrt(1 + ^(u[2],2))*(-u[5] + u[3] + u[4]))) - (dtdmp*u[6] + dtdtp*u[1])*(R*tauDiff*^(u[2],2)*sqrt(1 + ^(u[2],2))*(2*tauB*visc - 3*tauS*zeta + tauB*tauS*(u[3] - 2*u[4])) + tau*(-3*tauS*^(u[2],2)*(tauDiff*u[2]*zeta - R*sqrt(1 + ^(u[2],2))*zeta + R*tauDiff*u[5]) + tauB*tauDiff*u[2]*(3*R*u[2]*(u[3] + u[4]) - 3*tauS*(2*u[3] + u[4]) + ^(u[2],2)*(2*visc + tauS*(-2*u[3] + u[4]))) + R*tauB*sqrt(1 + ^(u[2],2))*(3*tauS*(-u[5] + u[3] + u[4]) + ^(u[2],2)*(4*visc - tauS*(3*u[5] + u[3] + u[4])))))))/(R*tau*sqrt(1 + ^(u[2],2))*(3*dtdtp*dtp*^(u[1],2)*tauB*tauS - 3*^(dtp,2)*u[1]*tauB*tauS*^(u[2],2) + 3*dtdtp*dtp*^(u[1],2)*tauB*tauS*^(u[2],2) + 3*dmp*u[6]*tauB*tauS*(dtdmp*u[6] + dtdtp*u[1] + (-dtp + dtdmp*u[6] + dtdtp*u[1])*^(u[2],2)) - 4*dtdtp*u[1]*tauB*^(u[2],2)*visc - 3*dtdtp*u[1]*tauS*^(u[2],2)*zeta + 3*dtdtp*u[1]*tauB*tauS*u[5] - 3*dtp*tauB*tauS*^(u[2],2)*u[5] + 3*dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[5] - 3*dtdtp*u[1]*tauB*tauS*u[3] + 3*dtp*tauB*tauS*^(u[2],2)*u[3] + dtdtp*u[1]*tauB*tauS*^(u[2],2)*u[3] + tauB*tauS*(3*dtp*^(u[2],2) + dtdtp*u[1]*(-3 + ^(u[2],2)))*u[4] + dtdmp*u[6]*(3*dtp*u[1]*tauB*tauS*(1 + ^(u[2],2)) - 3*tauS*^(u[2],2)*zeta + tauB*(3*tauS*u[5] - 3*tauS*(u[3] + u[4]) + ^(u[2],2)*(-4*visc + tauS*(3*u[5] + u[3] + u[4]))))))
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
    =#
    