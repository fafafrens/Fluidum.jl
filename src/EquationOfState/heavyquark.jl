
#struttura con due tipi di variabili, T e S
Base.@kwdef struct Heavy_Quark{T,S} <:EquationOfState 
    mass::T = 1.5 #campo (field) massa, di tipo T
    a1::T  = -15.526548963383643 #campo (field) a1, di tipo T
    a2::T  = 18.6159584620131
    b1::T  = -3.3147904483107595
    b2::T  = 5.310983721567554 
    a3::T  = -10.731808109698516
    a4::T  = 2.7413302179949346
    b4::T  = 1.8600649533271152
    b3::T  = -4.653922019495976 
    c::T  = -1.0465330501411811
    d::T  = 0.09551531822245873
    hadron_list::HadronResonaceGas{S}=HadronResonaceGas() #campo hadron_list, di tipo HadronResonanceGas{S}, 
    #a cui viene assegnato come valore HadronResonaceGas() (vedi charm_HRG.jl)
    #hadron_list field is an instance (istanza) of HadronResonaceGas with a specific type S
    #particle_list::HadronResonaceGas{S}=HadronResonaceGas() #campo hadron_list, di tipo HadronResonanceGas, 
    
end


@inline @fastmath function pressure(T,x::Heavy_Quark) 
    (exp((-0.01*^(x.d,2) - 0.1*^(x.c,2)*T)/^(T,2))*(0.0005624486560000001*x.a4 
    + T*(0.003652264*x.a3 + T*(0.023716000000000004*x.a2 + 0.154*x.a1*T + 
    5.208957878352717*^(T,2)))))/(0.0005624486560000001*x.b4 + 
    0.003652264*x.b3*T + 0.023716*x.b2*^(T,2) + 0.154*x.b1*^(T,3) + ^(T,4))*^(T,4)*fmGeV^3
end #T^4 = GeV^4

#@inline pressure(T,x::FluidProperties)=pressure(T,x.equation_of_state)

@inline function pressure_derivative(T,::Val{1},x::Heavy_Quark) 

    (-1.2653939625448257e-6*exp((-0.01*^(x.d,2) - 0.1*^(x.c,2)*T)/^(T,2))*(-633936.5739133607*^(T,2)*(0.014609056000000004*x.b4 + T*(0.071148*x.b3 + T*(0.30800000000000005*x.b2 + 
    1*x.b1*T)))*(0.0001079771941211913*x.a4 + T*(0.0007011506111765668*x.a3 + T*(0.004552926046601084*x.a2 + 0.02956445484805898*x.a1*T + 1*^(T,2)))) 
    + 121701.22867529064*^(T,2)*(0.014609056000000004*x.a4 + T*(0.071148*x.a3+ T*(0.30800000000000005*x.a2 + 1*x.a1*T)))*(0.0005624486560000001*x.b4 + 
    T*(0.003652264*x.b3 + T*(0.023716000000000004*x.b2 + 0.154*x.b1*T + 1*^(T,2)))) - 1.6465885036710668e7*^(T,2)*(0.0001079771941211913*x.a4 
    + T*(0.0007011506111765668*x.a3 + T*(0.004552926046601084*x.a2 + 0.02956445484805898*x.a1*T + 1*^(T,2))))*(0.0005624486560000001*x.b4 + 
    T*(0.003652264*x.b3 + T*(0.023716000000000004*x.b2 + 0.154*x.b1*T +1*^(T,2)))) - 82329.42518355335*(1*^(x.d,2) + 
    5*^(x.c,2)*T)*(0.0001079771941211913*x.a4 + T*(0.0007011506111765668*x.a3+ T*(0.004552926046601084*x.a2 + 0.02956445484805898*x.a1*T + 
    1*^(T,2))))*(0.0005624486560000001*x.b4 + T*(0.003652264*x.b3 + T*(0.023716000000000004*x.b2 + 0.154*x.b1*T + 1*^(T,2))))))/(^(1 + 
    (0.0005624486560000001*x.b4)/^(T,4) + (0.003652264*x.b3)/^(T,3) + (0.023716*x.b2)/^(T,2) + (0.154*x.b1)/T,2)*^(T,7))*fmGeV^3

end #T^3 = GeV^3

#pressure_derivative(T,::Val{1},x::FluidProperties )=pressure_derivative(T,Val(1),x.equation_of_state)

@inline function pressure_derivative(T,::Val{2},x::Heavy_Quark) 
    
    (exp((-0.01*^(x.d,2) - 0.1*^(x.c,2)*T)/^(T,2))*(-32*(2*x.b4 +     T*(11.36363636363636*x.b3 + T*(63.24843987181649*x.b2 + 
    342.2534625098294*x.b1*T + 1777.9400649861268*^(T,2))))   *(1*x.b4 +     T*(6.493506493506493*x.b3 + T*(42.165626581210994*x.b2 + 
    273.80277000786356*x.b1*T + 1777.940064986127*^(T,2))))*(^(T,2)*(1*x.a4    + T*(4.870129870129869*x.a3 + T*(21.082813290605497*x.a2 + 
    68.45069250196589*x.a1*T))) - 46.30607454374214*(1*^(x.d,2) +     5*^(x.c,2)*T)*(0.0001079771941211913*x.a4 + T*(0.0007011506111765668*x.a3 
    + T*(0.004552926046601084*x.a2 + 0.02956445484805898*x.a1*T +     1*^(T,2))))) +     3.5130428054668146e11*^(T,2)*(0.0001079771941211913*x.a4 +     T*(0.0007011506111765668*x.a3 + T*(0.004552926046601084*x.a2 + 
    0.02956445484805898*x.a1*T + 1*^(T,2))))*(1.476292956302297e-6*^(x.b4,2)     + x.b4*T*(0.000016433687825257478*x.b3 + T*(0.00008670371011702403*x.b2 + 
    0.0004186492829493334*x.b1*T + 0.0016873459680000003*^(T,2))) +     ^(T,2)*(0.000046686613139936006*^(x.b3,2) +     x.b3*T*(0.0005052663759733334*x.b2 + (0.0025310189520000004*x.b1 +     0.010956792000000002*T)*T) + ^(T,2)*(0.0014061216400000005*^(x.b2,2) + 
    x.b2*(0.014609056000000002*x.b1 + 0.06719533333333334*T)*T +     ^(T,2)*(0.03952666666666667*^(x.b1,2) + 0.385*x.b1*T + 1*^(T,2))))) +     (20*^(1*x.b4 + T*(6.493506493506493*x.b3 + T*(42.165626581210994*x.b2 +     273.80277000786356*x.b1*T +     1777.940064986127*^(T,2))),2)*(-2.7380277000786357*^(T,2)*(0.2*^(x.d,2) 
    + 1*^(x.c,2)*T)*(0.014609056000000004*x.a4 + T*(0.071148*x.a3 +     T*(0.30800000000000005*x.a2 + 1*x.a1*T))) + ^(T,4)*(1*x.a4 +     T*(3.8961038961038956*x.a3 + T*(12.649687974363298*x.a2 +     27.380277000786357*x.a1*T))) - 2.815071769561368*(^(T,2)*(1*^(x.d,2) + 
    3.3333333333333335*^(x.c,2)*T) - 0.006666666666666667*^(1*^(x.d,2) +     5*^(x.c,2)*T,2))*(0.0010656921903157896*x.a4 +     T*(0.006920079157894737*x.a3 + T*(0.044935578947368424*x.a2 +     0.29178947368421054*x.a1*T +     9.869604401089358*^(T,2))))))/^(T,2)))/^(1*x.b4 +     6.493506493506493*x.b3*T + 42.165626581210994*x.b2*^(T,2) + 
    273.80277000786356*x.b1*^(T,3) + 1777.940064986127*^(T,4),3)*fmGeV^3

end #T^2 = GeV^2

function free_charm(T,fug,x::Heavy_Quark)
    m = x.mass
    b2 = besselkx(2,m/T)
    b1 = besselk1x(m/T)
    #b2 = besselk(2,m/T)
    #b1 = besselk(1,m/T)
    b3 = b1+4/(m/T)*b2

    ex = exp(fug-m/T)
    deg = 6 #spin x2,color x3
    density = deg*(m^2* T /(2 *π^2)* ex* b2)* fmGeV^3; #fm-3
    densityDerT = deg*((ex*m^2*(m*b1 +2*T*b2 + m*b3))/(4*π^2*T))* fmGeV^3 ; #(*fm^-3/GeV*)
    densityDerμ= deg*(m^2* T /(2 *π^2)* ex* b2)* fmGeV^3+0.0001; #(*fm^-3/GeV*)
    densityDerTDerμ =deg*(ex*m^2*(m*b1 + 2*T*b2 + m*b3)/(4*π^2*T))* fmGeV^3 ;
    #@show b2, m, T
    return (density,densityDerT,densityDerμ,densityDerTDerμ )
end



struct HQdiffusion{M}<:Diffusion
    DsT::M
    mass::M
end

function HQdiffusion()
    HQdiffusion(0.39,1.5)  
end

struct QGPViscosity{T}<:ShearViscosity
    ηs::T
    Cs::T
end



@inline function bulk_viscosity(T,y::ZeroBulkViscosity)  
    zero(promote_type(typeof(T)))
end

@inline function τ_bulk(T,y::ZeroBulkViscosity)  
    one(promote_type(typeof(T)))
end


@inline function bulk_viscosity(T,entropy,y::SimpleBulkViscosity{N}) where {N}
    y.ζs/(1+((T-0.175)/0.024)^2)*invfmGeV*entropy
 end
 @inline function τ_bulk(T,entropy,dtdtp,y::SimpleBulkViscosity{N}) where {N}
     cs2= entropy/(T*dtdtp)
     bulk_viscosity(T,entropy,y)/(T*entropy*y.Cζ)*1/(1/3-cs2)^2+0.1
 end
 

@inline function viscosity(T,entropy,y::QGPViscosity{N}) where {N}
    #entropy GeV^3
    y.ηs*entropy*invfmGeV
    #zero(promote_type(typeof(T),typeof(entropy)))
 end


 
@inline function τ_shear(T,entropy,y::QGPViscosity{N}) where {N}
    viscosity(T,entropy,y)/((T*entropy)*y.Cs)
    #one(promote_type(typeof(T),typeof(entropy)))
end

function diffusion(T,n,x::HQdiffusion{M}) where {M}
    density = n #fm-^3
    #κ = x.DsT/T*density/fmGeV #fm^-2
    κ = x.DsT/T/fmGeV
end

function τ_diffusion(T,x::HQdiffusion{M}) where {M}
    m = x.mass
    #@show m/T
    #@show T
  
    b2 = besselk(2,m/T)
    b1 = besselk(1,m/T)
    b3 = b1+4/(m/T)*b2
    b4 = b2 + 6/(m/T)*b3
    b5 = b3+8/(m/T)*b4
    
    #b3 = besselk(3,m/T)
    #b4 = besselk(4,m/T)
    #b5 = besselk(5,m/T)
   
    #[1/GeV]*fm GeV/0.2 = [fm]/fmGeV
  return ((2*π*x.DsT) *m^3/T^3/(96*π*T)*(2*b1 - 3*b3 +b5)/b2)/fmGeV; #
  
end


function τ_diffusion(T,x::ZeroDiffusion) 
    return 1.0
end

function diffusion(T,n,x::ZeroDiffusion) 
    return 0.0
end
