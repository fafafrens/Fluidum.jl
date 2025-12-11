

struct HadronResonaceGas{T,Bool} <:EquationOfState
    particle_list::T
    relativistic::Bool
end



struct Particle{A,B,C} #for particle.data file
    Name::A
    Mass::B
    Gamma::B
    Spin::B
    Isospin::B
    I3::B
    Nq::B
    Ns::B
    Naq::B
    Nas::B
    Nc::B
    Nac::B
    MC::C
end



#read in resonaces
function HadronResonaceGas(;name_file=root_particle_lists*"/particles.data",Maxmass=2.1,Minmass=0.05,condition=x->true,relativistic = false)
    data     =readdlm(name_file,comment_char='#',comments=true)
    names    =convert.(String,data[:,1])
    mass     =convert.(Float64,data[:,2])
    gamma    =convert.(Float64,data[:,3])
    spin     =convert.(Float64,data[:,4])
    iso_spin =convert.(Float64,data[:,5])
    I3       =convert.(Float64,data[:,6])
    Nq       =convert.(Float64,data[:,7])
    Ns       =convert.(Float64,data[:,8])
    Naq      =convert.(Float64,data[:,9])
    Nas      =convert.(Float64,data[:,10])
    Nc       =convert.(Float64,data[:,11])
    Nac      =convert.(Float64,data[:,12])
    MC       =convert.(Int64,data[:,13])

   
    fulllist=StructArray(Particle.(
        names    ,
        mass     ,
        gamma    ,
        spin     ,
        iso_spin ,
        I3       ,
        Nq       ,
        Ns       ,
        Naq      ,
        Nas      ,
        Nc       ,
        Nac      ,
        MC       
       ))
    filterlist=filter(x->(x.Mass<Maxmass&&x.Mass>Minmass&&x.Name != "de2000plb"&&x.Name !="de2000plu"&&condition(x)),fulllist)
   HadronResonaceGas(filterlist,relativistic)

end



#pretty printing
function Base.show(io::IO, z::HadronResonaceGas)
    min ,max = extrema(z.particle_list.Mass)
    print(io,"Hadron Resonace gas: ",length(z.particle_list)," particles with mass ⊆ ",min,"..",max," GeV" )
end


function Base.show(io::IO, ::MIME"text/plain", z::HadronResonaceGas) 
    min ,max = extrema(z.particle_list.Mass)
    print(io,"Hadron Resonace gas: ",length(z.particle_list)," particles with mass ⊆ ",min,"..",max,"GeV\n " )
    for part in z.particle_list
        print(io,part.Name," mass=",part.Mass,"\n")
    end 
end

#condition to exclude particles, which are already included in Walecka model (also include neutrons?)
#waleckacondition(x)=x.Name != "f00600zer"&&x.Name !="om0782zer"&& x.Name != "pr0938plu"  && x.Name !="pr0938plb" && x.Name !="ne0939zer" && x.Name != "ne0939zrb"#remove neutrons


@inline Base.getindex(elm::HadronResonaceGas, i::Int)= elm.particle_list[i]
@inline Base.eachindex(elm::HadronResonaceGas)=Base.eachindex(elm.particle_list)


"""
    Pressure and its derivatives of a Hadron Gas, for either relativistic (Bose-Einstein/Fermi-Dirac) or non-relativistic (Maxwell-Boltzmann) statistics,
    for non zero chemical potential 
"""
function thermodynamic(T,μ,x::HadronResonaceGas{L,B})  where {L,B}
    relativistic = x.relativistic
    if isless(T,zero(T))
        return Thermodynamic(zero(T),(zero(T),zero(T)),(zero(T),zero(T),zero(T)))
    end

    pressure=zero(T)
    p10=zero(T)
    p01=zero(T)
    
    p20=zero(T)
    p02=zero(T)
    p11=zero(T)

    for i in eachindex(x) 
        #loop over all particles 
        QB = round((x[i].Nq+x[i].Ns+x[i].Nc-x[i].Naq-x[i].Nas-x[i].Nac)/3.0)
        m=x[i].Mass
        reducemass=m/T
        Spin=x[i].Spin
        degeneracy=(2*Spin+1)

        if reducemass<500*one(T)
            if relativistic == true
                if iseven(Spin*2)
            
                    for n in 1:4 #Loop over mesons: Bose-Einstein statistics 
                        
                        b0n=besselk0(n*reducemass)
                        b1n=besselk1(n*reducemass)
                        b2n=b0n+2/(n*reducemass)*b1n
                        expn=exp(n* μ* QB/T )
                
                        pressure+=degeneracy*T^4/(2 *pi^2)* (reducemass)^2/n^2*expn*b2n 
                        p01+=degeneracy/(2 *pi^2*n)*expn*m^2*QB*T*b2n 
                        p10+= degeneracy/(2 *pi^2*n^2)*expn* m^2 * (b2n *(4*T-μ*QB*n)+ b1n*(m*n))
                        p20+= degeneracy/(2 *pi^2 *n^3 *T^2)*expn*m*(
                        m* n*((m^2+(μ* QB)^2)*n^2-6*n*T*(μ* QB)+12*T^2)*b0n
                        +(m^2*n^2*(-2*(μ* QB)*n+5*T)+2*T*((μ* QB)^2*n^2-6*(μ* QB)*n*T+12*T^2 ))*b1n
                        )   
                        p11+= -degeneracy/(2 *pi^2*n*T)*expn*m^2*QB*(
                        -m*n*b1n +(μ* QB*n-3*T)*b2n)
                        p02+= degeneracy/(2 *pi^2)*expn*m^2*QB^2*b2n
                    end    
                
                else
                     
                    for n in 1:3 #Loop over baryons: Fermi-Dirac statistics
                        b0n=besselk0(n*reducemass)
                        b1n=besselk1(n*reducemass)
                        b2n=b0n+2/(n*reducemass)*b1n
                        expn=exp(n* μ* QB/T )
                
                        pressure+= degeneracy*T^4 /(2 *pi^2)* (reducemass)^2/n^2*(-1)^(n+1) *expn *b2n 
                        p01+= (-1)^(n+1)*degeneracy/(2 *pi^2*n)*expn *b2n *m^2*QB*T
                        p10+= (-1)^(n+1)*degeneracy/(2 *pi^2*n^2)*expn * m^2 *(b2n *(4*T-μ*QB*n)+ 
                            b1n*(m*n))
                        p20+=(-1)^(n)*degeneracy/(2 *pi^2 *n^3 *T^2)*expn*m*(
                            -m* n*((m^2+(μ* QB)^2)*n^2-6*n*T*(μ* QB)+12*T^2)*b0n
                            -(m^2*n^2*(-2*(μ* QB)*n+5*T)+2*T*((μ* QB)^2*n^2-6*(μ* QB)*n*T+12*T^2 ))*b1n
                        )
                        p11+= (-1)^(n)*degeneracy/(2 *pi^2*n*T)*expn*m^2*QB*(
                        -m*n* b1n +(μ* QB*n-3*T)*b2n)
                        p02+= (-1)^(n+1)*degeneracy/(2 *pi^2)*m^2*QB^2*expn *b2n 
                    end    
                end
            
            else #Maxwell-Boltzmann statistics 
                b0 = besselkx(0,reducemass)         
                b1 = besselk1x(reducemass)
                b2 = besselkx(2,reducemass)
                b3 = b1+4/(reducemass)*b2
                b4 = b2+6/(reducemass)*b3

                factor1 = m^2/2*(b0 + b4)
                factor2 = 2*m*(T-μ)*(b1 + b3)
                factor3 = (-4*μ*T+2*μ^2+m^2+4*T^2)*b2

                ex=exp(QB* μ - reducemass)
                
                pressure += degeneracy*(T^2*m^2 /(2 *π^2)* ex* b2) 
                
                p10 += degeneracy*((ex*m^2*(m*b1 +2*(2*T-μ)*b2 + m*b3))/(4*π^2)) #dP/dT
                p01 += degeneracy*(T*m^2 /(2 *π^2)* ex* b2) #dP/dmu
                p20 += degeneracy*((ex*m^2/T^2*(factor1+factor2+factor3))/(4*π^2)) #d2P/dT^2
                p02 += degeneracy*(m^2 /(2 *π^2)* ex* b2) #d2P/dmu^2
                p11 += degeneracy*((ex*m^2/T*(m*b1 +2*(T-μ)*b2 + m*b3))/(4*π^2)) #d2P/dmudT

            end
        else

            pressure+= zero(T)
            p01+= zero(T)
            p10+= zero(T)
            p20+=zero(T)
            p11+= zero(T)
            p02+= zero(T)

        end
    end

    return Thermodynamic(pressure*fmGeV3,(p10*fmGeV3,p01*fmGeV3),(p20*fmGeV^2,p11*fmGeV^2,p02*fmGeV^2))
end

"""
    Pressure and its derivatives of a Hadron Gas, for either relativistic (Bose-Einstein/Fermi-Dirac) or non-relativistic (Maxwell-Boltzmann) statistics, 
    for zero chemical potential 
"""
function thermodynamic(T,x::HadronResonaceGas{L,B})  where {L,B}
    μ = zero(T)
    return  thermodynamic(T,μ,x)
end