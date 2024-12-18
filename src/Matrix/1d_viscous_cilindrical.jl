# the convention here are T, ux,  \[Pi]yy, \[Pi]zz, \[Pi]B




function matrxi1d_visc!(A_i,Source,ϕ,t,X,params)
 
    p = pressure(ϕ[1],params.eos) #entropy
    dpt = pressure_derivative(ϕ[1],Val(1),params.eos) #entropy
    dptt = pressure_derivative(ϕ[1],Val(2),params.eos)
        
    etaVisc=viscosity(ϕ[1],dpt,params.shear)
    tauS=τ_shear(ϕ[1],dpt,params.shear)
    tauB=τ_bulk(ϕ[1],dpt,dptt,params.bulk)
    zeta=bulk_viscosity(ϕ[1],dpt,params.bulk)
    
 
    (At,Ax, source)=one_d_viscous_matrix(ϕ,t,X[1],p,dpt,dptt,zeta,etaVisc,tauS,tauB)
    

        
        Ainv= inv(At)

    
        A_mul_B!(A_i[1], Ainv,Ax)
        
        
        jgemvavx!(Source, Ainv,source)
    
        
    
    end 



function matrxi1d_visc(ϕ,t,X,params)
 

    p = pressure(ϕ[1],Val(1),params.eos) #entropy
    dpt = pressure_derivative(ϕ[1],Val(1),params.eos) #entropy
    dptt = pressure_derivative(ϕ[1],Val(2),params.eos)
            
    etaVisc=viscosity(ϕ[1],dpt,params.shear)
    tauS=τ_shear(ϕ[1],dpt,params.shear)
    tauB=τ_bulk(ϕ[1],dpt,dptt,params.bulk)
    zeta=bulk_viscosity(ϕ[1],dpt,params.bulk)
        
     
        (At,Ax, source)=one_d_viscous_matrix(ϕ,t,X[1],p,dpt,dptt,zeta,etaVisc,tauS,tauB)
        
    
            
        Ainv= inv(At)
    
        
          
        return (Ainv*source,Ainv*Ax)
            
        
end 



@inbounds function one_d_viscous_matrix(u,tau,R,p,dtp,dtdtp,zeta,visc,tauS,tauB)
    
At=SMatrix{5,5}(

    dtdtp*u[1] + (dtp + dtdtp*u[1])*^(u[2],2)

    ,(dtp + dtdtp*u[1])*u[2]*sqrt(1 + ^(u[2],2))

    ,0

    ,0

    ,0

    ,2*u[2]*(dtp*u[1] + u[5] - u[3] - u[4])

    ,((1 + 2*^(u[2],2))*(dtp*u[1] + u[5] - u[3] - u[4]))/sqrt(1 + ^(u[2],2))

    ,(-2*u[2]*visc)/(3*^(R,2)*sqrt(1 + ^(u[2],2)))

    ,(-2*u[2]*visc)/(3*^(tau,2)*sqrt(1 + ^(u[2],2)))

    ,(u[2]*zeta)/sqrt(1 + ^(u[2],2))

    ,-^(u[2],2)

    ,-(u[2]*sqrt(1 + ^(u[2],2)))

    ,(tauS*sqrt(1 + ^(u[2],2)))/^(R,2)

    ,0

    ,0

    ,-^(u[2],2)

    ,-(u[2]*sqrt(1 + ^(u[2],2)))

    ,0

    ,(tauS*sqrt(1 + ^(u[2],2)))/^(tau,2)

    ,0

    ,^(u[2],2)

    ,u[2]*sqrt(1 + ^(u[2],2))

    ,0

    ,0

    ,tauB*sqrt(1 + ^(u[2],2))

)
    #########################################################################################################################################################################
    
Ax=SMatrix{5,5}(

    (dtp + dtdtp*u[1])*u[2]*sqrt(1 + ^(u[2],2))

    ,dtp + (dtp + dtdtp*u[1])*^(u[2],2)

    ,0

    ,0

    ,0

    ,((1 + 2*^(u[2],2))*(dtp*u[1] + u[5] - u[3] - u[4]))/sqrt(1 + ^(u[2],2))

    ,2*u[2]*(dtp*u[1] + u[5] - u[3] - u[4])

    ,(-2*visc)/(3*^(R,2))

    ,(-2*visc)/(3*^(tau,2))

    ,zeta

    ,-(u[2]*sqrt(1 + ^(u[2],2)))

    ,-1 - ^(u[2],2)

    ,(tauS*u[2])/^(R,2)

    ,0

    ,0

    ,-(u[2]*sqrt(1 + ^(u[2],2)))

    ,-1 - ^(u[2],2)

    ,0

    ,(tauS*u[2])/^(tau,2)

    ,0

    ,u[2]*sqrt(1 + ^(u[2],2))

    ,1 + ^(u[2],2)

    ,0

    ,0

    ,tauB*u[2]

)
    #########################################################################################################################################################################

source=SVector{5}(

    (dtp*u[1]*(R + R*^(u[2],2) + tau*u[2]*sqrt(1 + ^(u[2],2))) + tau*u[2]*sqrt(1 + ^(u[2],2))*(u[5] - u[3] - u[4]) + R*((1 + ^(u[2],2))*u[5] + u[4] - ^(u[2],2)*(u[3] + u[4])))/(R*tau)

    ,(dtp*u[1]*u[2]*(tau*u[2] + R*sqrt(1 + ^(u[2],2))) + R*u[2]*sqrt(1 + ^(u[2],2))*(u[5] - u[3] - u[4]) + tau*(-2*u[3] + ^(u[2],2)*(u[5] - u[3] - u[4]) - u[4]))/(R*tau)

    ,(4*tau*u[2]*visc - 2*R*sqrt(1 + ^(u[2],2))*visc + 3*R*tau*u[3])/(3*^(R,3)*tau)

    ,(-2*tau*u[2]*visc + 4*R*sqrt(1 + ^(u[2],2))*visc + 3*R*tau*u[4])/(3*R*^(tau,3))

    ,(u[2]*zeta)/R + (sqrt(1 + ^(u[2],2))*zeta)/tau + u[5]

)

return (At,Ax, source)
end


