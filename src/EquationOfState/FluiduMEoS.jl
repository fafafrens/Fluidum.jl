
struct FluiduMEoS{T} <:EquationOfState
    a1::T
    a2::T
    a3::T
    a4::T
    b1::T
    b2::T
    b3::T
    b4::T
    c::T
    d::T

end

function FluiduMEoS()

    a1 = -15.526548963383643
    a2 = 18.6159584620131
    b1 = -3.3147904483107595
    b2 = 5.310983721567554 
    a3 = -10.731808109698516
    a4 = 2.7413302179949346
    b4 = 1.8600649533271152
    b3 = -4.653922019495976 
    c = -1.0465330501411811
    d = 0.09551531822245873
    FluiduMEoS(a1,a2,a3,a4,b1,b2,b3,b4,c,d)
end




@inline @fastmath function den(t,b1,b2,b3,b4)
    invt=1/t
    x=154*invt
    
    invt2=invt*invt
    invt3=invt*invt2
    fun=@evalpoly(x,1,b1,b2,b3,b4)
    fun1=@evalpoly(x,b1,2*b2,3*b3,4*b4)
    fun2=@evalpoly(x,2*b2,6*b3,12*b4)
    
    (fun,@muladd(-154*invt2*fun1) ,@muladd(308*invt3*fun1+(23716)*fun2*invt3*invt))
end

@inline @fastmath function den_pert(t,b1,b2,b3,b4)
    invt=1/t
    x=154*invt
    
    invt2=invt*invt
    invt3=invt*invt2
    invt4=invt*invt3
    fun=@evalpoly(x,1,b1,b2,b3,b4)
    fun1=@evalpoly(x,b1,2*b2,3*b3,4*b4)
    fun2=@evalpoly(x,2*b2,6*b3,12*b4)
    fun3=@evalpoly(x,6*b3,24*b4)
    
    (fun,@muladd(-154*invt2*fun1) ,@muladd(308*invt3*fun1+(23716)*fun2*invt3*invt),
    @muladd(-924*invt4*fun1-142_296*invt4*invt*fun2-3_652_264*invt4*invt2*fun3 )
    )
end

@inline @fastmath function num(t,a1,a2,a3,a4)
    
    invt=1/t
    x=154*invt
    invt2=invt*invt
    invt3=invt*invt2
    # ((95*pi^2)/180)=5.208957878352717
    fun=@evalpoly(x,5.208957878352717,a1,a2,a3,a4)
    fun1=@evalpoly(x,a1,2*a2,3*a3,4*a4)
    fun2=@evalpoly(x,2*a2,6*a3,12*a4)
    (fun,@muladd(-154*invt2*fun1) ,@muladd(308*invt3*fun1+(23716)*fun2*invt3*invt))
end

@inline @fastmath function num_pert(t,a1,a2,a3,a4)
    
    invt=1/t
    x=154*invt
    invt2=invt*invt
    invt3=invt*invt2
    invt4=invt*invt3
    # ((95*pi^2)/180)=5.208957878352717
    fun=@evalpoly(x,5.208957878352717,a1,a2,a3,a4)
    fun1=@evalpoly(x,a1,2*a2,3*a3,4*a4)
    fun2=@evalpoly(x,2*a2,6*a3,12*a4)
    fun3=@evalpoly(x,6*a3,24*a4)
    
    (fun,@muladd(-154*invt2*fun1) ,@muladd(308*invt3*fun1+(23716)*fun2*invt3*invt),
    @muladd(-924*invt4*fun1-142_296*invt4*invt*fun2-3_652_264*invt4*invt2*fun3 )
    )
end

#Exp[-c^2/(t/100)-d^2/((t)/100)^2]


@inline @fastmath function exponent(t,c,d)
    
    c1=-c^2 
    d1=-d^2 
    invt=1/t
    x=100*invt
    invt2=invt*invt
    invt3=invt*invt2
    fac1=exp(muladd(d1,x,c1)*x)
    #fac1=exp(c1*x+d1*(x)*x) 
    tinypol=muladd(x,2*d1,c1)
    
    der=fac1*tinypol
    der2=muladd(tinypol,der ,fac1*(2*d1))



    (fac1,@muladd(-100*invt2*der) ,@muladd(200*invt3*der+(10_000)*der2*invt3*invt))
    #(fac1,der,der2)

end

@inline @fastmath function exponent_pert(t,c,d)
    
    c1=-c^2 
    d1=-d^2 
    invt=1/t
    x=100*invt
    invt2=invt*invt
    invt3=invt*invt2
    invt4=invt3*invt
    fac1=exp(muladd(d1,x,c1)*x)
    #fac1=exp(c1*x+d1*(x)*x) 
    tinypol=muladd(x,2*d1,c1)
    
    der=fac1*tinypol
    der2=muladd(tinypol,der ,fac1*(2*d1))
    @muladd der3=6*d1*der+tinypol*tinypol*der
    #fac1*((c1+2*d1*x)^2 +(2*d1))


    (fac1,@muladd(-100*invt2*der) ,@muladd(200*invt3*der+(10_000)*der2*invt3*invt),
    @muladd(-600*invt4*der-60_000*invt4*invt*der2-1_000_000*invt4*invt2*der3))
    #(fac1,der,der2)

end

@inline @fastmath function thermodynamic(T,x::FluiduMEoS)
    
    t=T*1000
    (d,d1,d2)=den(t,x.b1,x.b2,x.b3,x.b4)
    (n,n1,n2)=num(t,x.a1,x.a2,x.a3,x.a4)
    (e,e1,e2)=exponent(t,x.c,x.d)
    invd=1/d
    invd2=invd*invd
    invd3=invd2*invd
    @muladd fun=e*n*invd*fmGeV3
    @muladd fun1=(d*n*e1+e*(d*n1-n*d1))*invd2*fmGeV3

    @muladd fun2=(d*(2*e1*(-n*d1+d*n1)+d*n*e2)+e*(n*(2*d1*d1-d*d2)+d*(-2*d1*n1+d*n2)))*invd3*fmGeV3
    #@muladd fun3= (d*((6*d1*d1*e1 - 3*d*d2*e1 - 3*d*d1*e2 + d*d*e3)*n + 3*d*(-2*d1*e1*n1 + d*e2*n1 + d*e1*n2)) + e*(-((6*d1*d1*d1 - 6*d*d1*d2 + d*d*d3)*n) + d*(6*d1*d1*n1 - 3*d*d2*n1 - 3*d*d1*n2 + d*d*n3)))*invd4
    T2=T*T
    T3=T2*T
    zero_d=T3*T*fun
    @muladd one_d=T3*(4*fun+1_000*T*fun1)
    @muladd two_d=T2*(12*fun+T*(8_000*fun1+1_000_000*T*fun2))

    #@mullad three_d= T*(24*fun +36*t*fun1 +12*t*t*fun2+t*t*t*fun3)
    
    Thermodynamic(zero_d,(one_d,),(two_d,))

end

@inline @fastmath function thermodynamic_perturbation(T,x::FluiduMEoS)
    
    t=T*1000
    (d,d1,d2,d3)=den_pert(t,x.b1,x.b2,x.b3,x.b4)
    (n,n1,n2,n3)=num_pert(t,x.a1,x.a2,x.a3,x.a4)
    (e,e1,e2,e3)=exponent_pert(t,x.c,x.d)
    invd=1/d
    invd2=invd*invd
    invd3=invd2*invd
    invd4=invd3*invd
    @muladd fun=e*n*invd*fmGeV3
    @muladd fun1=(d*n*e1+e*(d*n1-n*d1))*invd2*fmGeV3

    @muladd fun2=(d*(2*e1*(-n*d1+d*n1)+d*n*e2)+e*(n*(2*d1*d1-d*d2)+d*(-2*d1*n1+d*n2)))*invd3*fmGeV3
    @muladd fun3= (d*((6*d1*d1*e1 - 3*d*d2*e1 - 3*d*d1*e2 + d*d*e3)*n + 3*d*(-2*d1*e1*n1 + d*e2*n1 + d*e1*n2)) + e*(-((6*d1*d1*d1 - 6*d*d1*d2 + d*d*d3)*n) + d*(6*d1*d1*n1 - 3*d*d2*n1 - 3*d*d1*n2 + d*d*n3)))*invd4*fmGeV3
    
    T2=T*T
    T3=T2*T
    zero_d=T3*T*fun
    @muladd one_d=T3*(4*fun+t*fun1)
    @muladd two_d=T2*(12*fun+t*(8*fun1+t*fun2))
    @muladd three_d=T*(24*fun+36*t*fun1+12*t*t*fun2+t*t*t*fun3)
    
    ThermodynamicPerturbation(zero_d,(one_d,),(two_d,),(three_d,))

end