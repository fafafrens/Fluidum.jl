# the convention here are T, ux,  \[Pi]yy, \[Pi]zz, \[Pi]B
# the shear rows carry the DNMR term dpp pi theta, dpp = delta_pipi(tauS) = 4/3 tauS by default


"""
    matrix1d_visc!(; delta_pipi = true)
1+1D viscous hydro; `delta_pipi = false` drops the DNMR term (plain Israel-Stewart).
"""
function matrix1d_visc!(; delta_pipi::Bool = true)
    (A_i, Source, ϕ, t, X, params) -> matrix1d_visc!(A_i, Source, ϕ, t, X, params; delta_pipi = delta_pipi)
end


function matrix1d_visc!(A_i,Source,ϕ,t,X,params; delta_pipi::Bool = true)
 
    p = pressure(ϕ[1],params.eos) #entropy
    dpt = pressure_derivative(ϕ[1],Val(1),params.eos) #entropy
    dptt = pressure_derivative(ϕ[1],Val(2),params.eos)
        
    etaVisc=viscosity(ϕ[1],dpt,params.shear)
    tauS=τ_shear(ϕ[1],dpt,params.shear)
    tauB=τ_bulk(ϕ[1],dpt,dptt,params.bulk)
    zeta=bulk_viscosity(ϕ[1],dpt,params.bulk)
    
 
    (At,Ax, source)= delta_pipi ? one_d_viscous_matrix(ϕ,t,X[1],p,dpt,dptt,zeta,etaVisc,tauS,tauB) :
                                  one_d_viscous_matrix(ϕ,t,X[1],p,dpt,dptt,zeta,etaVisc,tauS,tauB,zero(tauS))
    

        
        Ainv= inv(At)

    
        A_mul_B!(A_i[1], Ainv,Ax)
        
        
        jgemvavx!(Source, Ainv,source)
    
        
    
    end 



function matrix1d_visc(ϕ,t,X,params)
 

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



@inbounds function one_d_viscous_matrix(u,tau,R,p,dtp,dtdtp,zeta,visc,tauS,tauB,dpp=delta_pipi(tauS))
    
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

iszero(dpp) && return (At,Ax, source)
    # dpp pi theta in the two shear rows
dAt=SMatrix{5,5}(
    0,0,0,0,0
    ,0,0,(u[3]*u[2])/(^(R,2)*sqrt(1 + ^(u[2],2))),(u[4]*u[2])/(^(tau,2)*sqrt(1 + ^(u[2],2))),0
    ,0,0,0,0,0
    ,0,0,0,0,0
    ,0,0,0,0,0
)
dAx=SMatrix{5,5}(
    0,0,0,0,0
    ,0,0,u[3]/^(R,2),u[4]/^(tau,2),0
    ,0,0,0,0,0
    ,0,0,0,0,0
    ,0,0,0,0,0
)
dsource=SVector{5}(
    0,0
    ,(u[3]*(sqrt(1 + ^(u[2],2))/tau + u[2]/R))/^(R,2)
    ,(u[4]*(sqrt(1 + ^(u[2],2))/tau + u[2]/R))/^(tau,2)
    ,0
)
return (At + dpp*dAt, Ax + dpp*dAx, source + dpp*dsource)
end

# the DNMR coefficient of the shear rows
@inline delta_pipi(tauS) = (4/3) * tauS

