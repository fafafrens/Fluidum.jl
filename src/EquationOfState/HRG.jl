

struct HadronResonaceGas{T} <:EquationOfState
    particle_list::T
end



struct HadronResonaceGasNew{T} <:EquationOfState
    particle_list::T
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

struct NewParticle{A,B,C} #for therminator type of file, e.g. PDG2016Plus_massorder.dat
        ID::A
        Name::B 
       Mass ::C
       Width ::C
       Degeneracy::A 
       Baryon::A
       Strangeness::A
       Charm::A
       Bottom::A
       Isospin::C
       ElectricCharge::A
       N_decay_channels::A
end




#read in resonaces
function HadronResonaceGas(;name_file=artifacts"particle_lists"*"/particles.data",Maxmass=2.1,Minmass=0.05,condition=x->true)
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
   HadronResonaceGas(filterlist)

end


#read in resonances
function HadronResonaceGasNew(;name_file=artifacts"particle_lists"*"/PDG2016Plus_massorder.dat",Maxmass=2.1,Minmass=0.05,condition=x->true)
    data     =readdlm(name_file,comment_char='#',comments=true)
    ID          =convert.(Int64,data[:,1])
    Name =convert.(String,data[:,2])
    Mass =convert.(Float64,data[:,3])
    Width =convert.(Float64,data[:,4])
    Degeneracy =convert.(Int64,data[:,5]) 
    Baryon =convert.(Int64,data[:,6])
    Strangeness=convert.(Int64,data[:,7])
    Charm=convert.(Int64,data[:,8])
    Bottom=convert.(Int64,data[:,9])
    Isospin=convert.(Float64,data[:,10])
    ElectricCharge=convert.(Int64,data[:,11])
    N_decay_channels=convert.(Int64,data[:,12])

    fulllist=StructArray(NewParticle.(
       ID,
       Name, 
       Mass ,
       Width ,
       Degeneracy, 
       Baryon,
       Strangeness,
       Charm,
       Bottom,
       Isospin,
       ElectricCharge,
       N_decay_channels,
       ))
    filterlist=filter(x->(x.Mass<Maxmass&&x.Mass>Minmass&&condition(x)),fulllist)
   HadronResonaceGasNew(filterlist)

end

#pretty printing
function Base.show(io::IO, z::HadronResonaceGas)
    min ,max = extrema(z.particle_list.Mass)
    print(io,"Hadron Resonace gas: ",length(z.particle_list)," particles with mass ⊆ ",min,"..",max," GeV" )
end


function Base.show(io::IO, ::MIME"text/plain", z::HadronResonaceGas) where{T}
    min ,max = extrema(z.particle_list.Mass)
    print(io,"Hadron Resonace gas: ",length(z.particle_list)," particles with mass ⊆ ",min,"..",max,"GeV\n " )
    for part in z.particle_list
        print(io,part.Name," mass=",part.Mass,"\n")
    end 
end

#condition to exclude particles, which are already included in Walecka model (also include neutrons?)
waleckacondition(x)=x.Name != "f00600zer"&&x.Name !="om0782zer"&& x.Name != "pr0938plu"  && x.Name !="pr0938plb" && x.Name !="ne0939zer" && x.Name != "ne0939zrb"#remove neutrons
#prova=HadronResonaceGas(name_file="EquationofState/particles.data",condition=waleckacondition)





@inline Base.getindex(elm::HadronResonaceGas, i::Int)= elm.particle_list[i]
@inline Base.eachindex(elm::HadronResonaceGas)=Base.eachindex(elm.particle_list)

@inline Base.getindex(elm::HadronResonaceGasNew, i::Int)= elm.particle_list[i]
@inline Base.eachindex(elm::HadronResonaceGasNew)=Base.eachindex(elm.particle_list)

"""
    thermodynamic(T,μ,x::HadronResonaceGas{L})  where {L}

TBW test test hello
"""
function thermodynamic(T,μ,x::HadronResonaceGas{L})  where {L}

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
        ## loop over all particles 
        QB = round((x[i].Nq+x[i].Ns+x[i].Nc-x[i].Naq-x[i].Nas-x[i].Nac)/3.0)
        #QS = round(x[i].Nas-x[i].Ns)
        #QC = round(x[i].Nc-x[i].Nac)
        m=x[i].Mass
        reducemass=m/T
        Spin=x[i].Spin
        
        if reducemass<500*one(T)
            if iseven(Spin*2)
          
                ## loop over mesons
                for n in 1:4
                    #besslk2n=besselk(2,n*reducemass) 
                    #@show T 
                    #@show m
                    #@show n*reducemass
                    besslk0n=besselk0(n*reducemass)
                    besslk1n=besselk1(n*reducemass)
                    #@show besslk0n
                    besslk2n=besslk0n+2/(n*reducemass)*besslk1n
                    expn=exp(n* μ* QB/T )
            
                    pressure+=T^4/(2 *pi^2)* (reducemass)^2*(2*Spin+1)/n^2*expn*besslk2n 
                    p01+=(2*Spin+1)/(2 *pi^2*n)*expn*m^2*QB*T*besslk2n 
                    p10+= (2*Spin+1)/(2 *pi^2*n^2)*expn* m^2 * (
                    besslk2n *(4*T-μ*QB*n)+ 
                    besslk1n*(m*n)
                    )
                    p20+= (2*Spin+1)/(2 *pi^2 *n^3 *T^2)*expn*m*(
                    m* n*((m^2+(μ* QB)^2)*n^2-6*n*T*(μ* QB)+12*T^2)*besslk0n
                    +(m^2*n^2*(-2*(μ* QB)*n+5*T)+2*T*((μ* QB)^2*n^2-6*(μ* QB)*n*T+12*T^2 ))*besslk1n
                    )   
                    p11+= -(2*Spin+1)/(2 *pi^2*n*T)*expn*m^2*QB*(
                    -m*n*besslk1n +(μ* QB*n-3*T)*besslk2n)
                    p02+= (2*Spin+1)/(2 *pi^2)*expn*m^2*QB^2*besslk2n
                end    
            else
                ## loop over barions 
                for n in 1:3
                    #@show T 
                    #@show m
                    #@show n*reducemass
                    besslk0n=besselk0(n*reducemass)
                    besslk1n=besselk1(n*reducemass)
                    besslk2n=besslk0n+2/(n*reducemass)*besslk1n
                    expn=exp(n* μ* QB/T )
            
                    pressure+= T^4 /(2 *pi^2)* (reducemass)^2* (2*Spin+1)/n^2*(-1)^(n+1) *expn *besslk2n 
                    p01+= (-1)^(n+1)*(2*Spin+1)/(2 *pi^2*n)*expn *besslk2n *m^2*QB*T
                    p10+= (-1)^(n+1)*(2*Spin+1)/(2 *pi^2*n^2)*expn * m^2 *(
                        besslk2n *(4*T-μ*QB*n)+ 
                        besslk1n*(m*n))
                    p20+=(-1)^(n)*(2*Spin+1)/(2 *pi^2 *n^3 *T^2)*expn*m*(
                        -m* n*((m^2+(μ* QB)^2)*n^2-6*n*T*(μ* QB)+12*T^2)*besslk0n
                        -(m^2*n^2*(-2*(μ* QB)*n+5*T)+2*T*((μ* QB)^2*n^2-6*(μ* QB)*n*T+12*T^2 ))*besslk1n
                    )
                    p11+= (-1)^(n)*(2*Spin+1)/(2 *pi^2*n*T)*expn*m^2*QB*(
                    -m*n* besslk1n +(μ* QB*n-3*T)*besslk2n)
                    p02+= (-1)^(n+1)*(2*Spin+1)/(2 *pi^2)*m^2*QB^2*expn *besslk2n 
                end
            
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

    return Thermodynamic(pressure,(p10,p01),(p20,p11,p02))
end

#prova=HadronResonaceGas(Maxmass=2.1)
#using BenchmarkTools


#@code_warntype thermodyanmics(0.1,0.1,prova)
#@benchmark thermodyanmics(0.1,0.1,$prova)

#@benchmark pressure(0.1,0.1,$prova)


"""
    thermodynamic(T,μ,x::HadronResonaceGasNew{L})  where {L}

TBW
"""
function thermodynamic(T,μ,x::HadronResonaceGasNew{L})  where {L}

    pressure=zero(T)
    p10=zero(T)
    p01=zero(T)
    
    p20=zero(T)
    p02=zero(T)
    p11=zero(T)

    for i in eachindex(x) 
        ## loop over all particles 
        QB = x[i].Baryon
        #QS = round(x[i].Nas-x[i].Ns)
        #QC = round(x[i].Nc-x[i].Nac)
        m=x[i].Mass
        reducemass=m/T
        degeneracy=x[i].Degeneracy
    
        if QB==0
          
            ## loop over mesons
            for n in 1:4
                besslk0n=besselk(0,n*reducemass)
                besslk1n=besselk(1,n*reducemass)
                besslk2n=besslk0n+2/(n*reducemass)*besslk1n
                expn=exp(n* μ* QB/T )
            
                pressure+=T^4/(2 *pi^2)* (reducemass)^2* degeneracy/n^2*expn*besslk2n 
                p01+=degeneracy/(2 *pi^2*n)*expn*m^2*QB*T*besslk2n 
                p10+= degeneracy/(2 *pi^2*n^2)*expn* m^2 * (
                    besslk2n *(4*T-μ*QB*n)+ 
                    besslk1n*(m*n)
                    )
                p20+= degeneracy/(2 *pi^2 *n^3 *T^2)*expn*m*(
                    m* n*((m^2+(μ* QB)^2)*n^2-6*n*T*(μ* QB)+12*T^2)*besslk0n
                    +(m^2*n^2*(-2*(μ* QB)*n+5*T)+2*T*((μ* QB)^2*n^2-6*(μ* QB)*n*T+12*T^2 ))*besslk1n
                    )   
                p11+= -degeneracy/(2 *pi^2*n*T)*expn*m^2*QB*(
                    -m*n*besslk1n +(μ* QB*n-3*T)*besslk2n)
                p02+= degeneracy/(2 *pi^2)*expn*m^2*QB^2*besslk2n
            end    
        else
            ## loop over barions 
            for n in 1:3
                besslk0n=besselk(0,n*reducemass)
                besslk1n=besselk(1,n*reducemass)
                besslk2n=besslk0n+2/(n*reducemass)*besslk1n
                expn=exp(n* μ* QB/T )
            
                pressure+= T^4 /(2 *pi^2)* (reducemass)^2* degeneracy/n^2*(-1)^(n+1) *expn *besslk2n 
                p01+= (-1)^(n+1)*degeneracy/(2 *pi^2*n)*expn *besslk2n *m^2*QB*T
                p10+= (-1)^(n+1)*degeneracy/(2 *pi^2*n^2)*expn * m^2 *(
                    besslk2n *(4*T-μ*QB*n)+ 
                    besslk1n*(m*n))
                p20+=(-1)^(n)*degeneracy/(2 *pi^2 *n^3 *T^2)*expn*m*(
                    -m* n*((m^2+(μ* QB)^2)*n^2-6*n*T*(μ* QB)+12*T^2)*besslk0n
                    -(m^2*n^2*(-2*(μ* QB)*n+5*T)+2*T*((μ* QB)^2*n^2-6*(μ* QB)*n*T+12*T^2 ))*besslk1n
                )
                p11+= (-1)^(n)*degeneracy/(2 *pi^2*n*T)*expn*m^2*QB*(
                    -m*n* besslk1n +(μ* QB*n-3*T)*besslk2n)
                p02+= (-1)^(n+1)*degeneracy/(2 *pi^2)*m^2*QB^2*expn *besslk2n 
            end
            
        end
    end
    return Thermodynamic(pressure,(p10,p01),(p20,p11,p02))
end
