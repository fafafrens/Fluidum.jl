struct HadronResonaceGas{T} <:EquationOfState
    particle_list::T
end

struct HadronResonaceGasNew{T} <:EquationOfState
    particle_list::T
end

struct HadronResonaceGas_ccbar{T,S} <:EquationOfState
    particle_list::T
    ccbar::S
end

struct Particle{A,B,C} #for particle.data file
    Name::A
    Mass::B
    Gamma::B
    Degeneracy::C
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




#read in resonances
function HadronResonaceGas(;name_file=string(root_particle_lists,"/OpenCharmParticleList_corrJS.txt"),Maxmass=4,Minmass=1.0,condition=x->true)
    data     =readdlm(name_file,comment_char='#',comments=true)
    names    =convert.(String,data[:,1])
    mass     =convert.(Float64,data[:,2])
    gamma    =convert.(Float64,data[:,3])
    deg      =convert.(Int64,data[:,4])
    spin     =convert.(Float64,data[:,5])
    iso_spin =convert.(Float64,data[:,6])
    I3       =convert.(Float64,data[:,7])
    Nq       =convert.(Float64,data[:,8])
    Ns       =convert.(Float64,data[:,9])
    Naq      =convert.(Float64,data[:,10])
    Nas      =convert.(Float64,data[:,11])
    Nc       =convert.(Float64,data[:,12])
    Nac      =convert.(Float64,data[:,13])
    MC       =convert.(Int64,data[:,14])

    fulllist=StructArray(Particle.(
        names    ,
        mass     ,
        gamma    ,
        deg      ,
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
    filterlist=filter(x->(x.Mass<Maxmass&&x.Mass>Minmass&&x.Name != "de2000plb"&&x.Name !="de2000plu" &&condition(x)),fulllist)
   HadronResonaceGas(filterlist)

end

function readresonancelist(;name_file=string(root_particle_lists,"/OpenCharmParticleList_corrJS.txt"),Maxmass=4,Minmass=1.0,condition=x->true)

    data     =readdlm(name_file,comment_char='#',comments=true)
    names    =convert.(String,data[:,1])
    mass     =convert.(Float64,data[:,2])
    gamma    =convert.(Float64,data[:,3])
    deg      =convert.(Int64,data[:,4])
    spin     =convert.(Float64,data[:,5])
    iso_spin =convert.(Float64,data[:,6])
    I3       =convert.(Float64,data[:,7])
    Nq       =convert.(Float64,data[:,8])
    Ns       =convert.(Float64,data[:,9])
    Naq      =convert.(Float64,data[:,10])
    Nas      =convert.(Float64,data[:,11])
    Nc       =convert.(Float64,data[:,12])
    Nac      =convert.(Float64,data[:,13])
    MC       =convert.(Int64,data[:,14])

    fulllist=StructArray(Particle.(
        names    ,
        mass     ,
        gamma    ,
        deg      ,
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
    return filter(x->(x.Mass<Maxmass&&x.Mass>Minmass&&x.Name != "de2000plb"&&x.Name !="de2000plu" &&condition(x)),fulllist)
   
end


#read in resonances
function HadronResonaceGasNew(;name_file=root_particle_lists*"/PDG2016Plus_massorder.dat",Maxmass=3.2,Minmass=1.8,condition=x->true)
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
    filterlist=filter(x->(x.Mass<Maxmass&&x.Mass>Minmass&&x.Charm!=0&&condition(x)),fulllist)
   HadronResonaceGasNew(filterlist)
end

#pretty printing
function Base.show(io::IO, z::HadronResonaceGas)
    min ,max = extrema(z.particle_list.Mass)
    print(io,"Hadron Resonace gas: ",length(z.particle_list)," particles with mass ⊆ ",min,"..",max," GeV" )
    for part in z.particle_list
        print(io,part.Name," mass=",part.Mass,"\t","spin=",part.Spin," QC=",(part.Nac+part.Nc)," deg=",(part.Degeneracy),"\n")
    end 
end


function Base.show(io::IO, ::MIME"text/plain", z::HadronResonaceGas) 
    min ,max = extrema(z.particle_list.Mass)
    print(io,"Hadron Resonace gas: ",length(z.particle_list)," particles with mass ⊆ ",min,"..",max,"GeV\n " )
    for part in z.particle_list
        print(io,part.Name," mass=",part.Mass,"\n")
    end 
end

function Base.show(io::IO, z::HadronResonaceGasNew)
    min ,max = extrema(z.particle_list.Mass)
    print(io,"Hadron Resonace gas: ",length(z.particle_list)," particles with mass ⊆ ",min,"..",max," GeV" )
    for part in z.particle_list
        print(io,part.Name," mass=",part.Mass,"\t","charm=",part.Charm,"\n")
    end 
end

function Base.show(io::IO, z::HadronResonaceGas_ccbar)
    min ,max = extrema(z.particle_list.Mass)
    print(io,"Hadron Resonace gas: ",length(z.particle_list)," particles with mass ⊆ ",min,"..",max," GeV" )
    for part in z.particle_list
        print(io,part.Name," mass=",part.Mass,"\n")
    end 
end


#condition to exclude particles, which are already included in Walecka model (also include neutrons?)
waleckacondition(x)=x.Name != "f00600zer"&&x.Name !="om0782zer"&& x.Name != "pr0938plu"  && x.Name !="pr0938plb" && x.Name !="ne0939zer" && x.Name != "ne0939zrb"#remove neutrons



@inline Base.getindex(elm::HadronResonaceGas, i::Int)= elm.particle_list[i] #to access the elements of particle_list
@inline Base.eachindex(elm::HadronResonaceGas)=Base.eachindex(elm.particle_list)

@inline Base.getindex(elm::HadronResonaceGasNew, i::Int)= elm.particle_list[i]
@inline Base.eachindex(elm::HadronResonaceGasNew)=Base.eachindex(elm.particle_list)

@inline Base.getindex(elm::HadronResonaceGas_ccbar, i::Int)= elm.particle_list[i]
@inline Base.eachindex(elm::HadronResonaceGas_ccbar)=Base.eachindex(elm.particle_list)


"""
Density of hadron resonance gas of charm 
"""
function thermodynamic(T,μ,x::HadronResonaceGas{L})  where {L}

    if isless(T,zero(T))
        return Thermodynamic(zero(T),(zero(T),zero(T)),(zero(T),zero(T),zero(T)))
    end

    density=zero(T)
    n10=zero(T)
    n01=zero(T)
    n20=zero(T)
    n02=zero(T)
    n11=zero(T)

    for i in eachindex(x) ## loop over all particles 
        QC = round((x[i].Nc+x[i].Nac))
        m=x[i].Mass
        
        reducemass=m/T
        degeneracy=x[i].Degeneracy
        
        if reducemass<500*one(T)
           
                b2 = besselkx(2,reducemass)
                b1 = besselk1x(reducemass)
                b3 = b1+4/(reducemass)*b2
                ex=exp(QC* μ - reducemass)
                fact = 1 #no correction factor in this case 
                
                density += fact*QC*degeneracy*(m^2* T /(2 *π^2)* ex* b2); #fm-3
                n10 += fact*QC*degeneracy*((ex*m^2*(m*b1 +2*T*b2 + m*b3))/(4*π^2*T)); #(*fm^-3/GeV*)
                n01+= fact*QC*QC*degeneracy*(m^2* T /(2 *π^2)* ex* b2); #(*fm^-3/GeV*)
                n11 +=fact*QC*QC*degeneracy*(ex*m^2*(m*b1 + 2*T*b2 + m*b3)/(4*π^2*T));
            
        else 
                density+= zero(T)
                n01+= zero(T)
                n10+= zero(T)
                n11+= zero(T)
        end
    end

    return Thermodynamic(density*fmGeV3,(n10*fmGeV3,n01*fmGeV3),(n20*fmGeV3,n11*fmGeV3,n02*fmGeV3))
end


"""
    Density of hadron resonance gas of charm, considering the canonical correction factor
"""
function thermodynamic(T,μ,x::HadronResonaceGas_ccbar{L,M})  where {L,M}
    ccbar = x.ccbar
    if isless(T,zero(T))
        return Thermodynamic(zero(T),(zero(T),zero(T)),(zero(T),zero(T),zero(T)))
    end

    density=zero(T)
    n10=zero(T)
    n01=zero(T)
    n20=zero(T)
    n02=zero(T)
    n11=zero(T)

    for i in eachindex(x) ## loop over all particles 
        QC = round((x[i].Nc+x[i].Nac))
        m=x[i].Mass
        
        reducemass=m/T
        degeneracy=x[i].Degeneracy
        
        if reducemass<500*one(T)
           
                b1 = besselk1x(reducemass)
                b2 = besselkx(2,reducemass)
                b3 = b1+4/(reducemass)*b2
                ex=exp(QC* μ - reducemass)
                
                #correction factor due to canonical ensemble 
                if QC == 1 
                    fact = besseli(1, ccbar/2)./besseli(0, ccbar/2)
                else 
                    fact = 1
                end   

                #density += fact*QC*degeneracy*(m^2* T /(2 *π^2)* ex* b2); #fm-3
                density += fact*QC*degeneracy*(ex*m^2*T*b2)/(2*pi^2); #fm-3



                
                n10 += fact*QC*degeneracy*((ex*m^2*(m*b1 +2*T*b2 + m*b3))/(4*π^2*T)); #dn/dT (fm^-3/GeV)
                n01+= QC*density; #dn/dfug (fm^-3/GeV)
                n11 +=fact*QC*QC*degeneracy*(ex*m^2*(m*b1 + 2*T*b2 + m*b3)/(4*π^2*T)); #dn/dTdfug (fm^-3/GeV)



                #n10 += fact*QC*degeneracy*((ex*m^2*(m*b1 +2*T*b2 + m*b3))/(4*π^2*T)); #(*fm^-3/GeV*)
                #n01+= fact*QC*QC*degeneracy*(m^2* T /(2 *π^2)* ex* b2); #(*fm^-3/GeV*)
                #n11 +=fact*QC*QC*degeneracy*(ex*m^2*(m*b1 + 2*T*b2 + m*b3)/(4*π^2*T));
        else 
                    density+= zero(T)
                    n01+= zero(T)
                    n10+= zero(T)
                    n11+= zero(T)

        end
    end

    return Thermodynamic(density*fmGeV3,(n10*fmGeV3,n01*fmGeV3),(n20*fmGeV3,n11*fmGeV3,n02*fmGeV3))
end

"""
    Density of hadron resonance gas of charm, considering relativistic (Bose-Einstein/Fermi-Dirac) dynamics  
"""
function thermodynamic(T,μ,x::HadronResonaceGasNew{L})  where {L}

    density=zero(T)
    n10=zero(T)
    n01=zero(T)
    n20=zero(T)
    n02=zero(T)
    n11=zero(T)

    for i in eachindex(x) 
        ## loop over all particles 
        QB = x[i].Baryon
        QC = round((x[i].Charm))
        m=x[i].Mass
        reducemass=m/T
        degeneracy=x[i].Degeneracy*2
        
    if reducemass<1000*one(T)
        if QB==0
          
            ## loop over mesons
            b2 = besselkx(2,reducemass)
            b1 = besselk1x(reducemass)
            b3 = b1+4/(reducemass)*b2
            ex=exp(QC* μ - reducemass)
            density += QC*degeneracy*(m^2* T /(2 *π^2)* ex* b2); #fm-3 added a QC term for J/Psi
            n10 += QC*degeneracy*((ex*m^2*(m*b1 +2*T*b2 + m*b3))/(4*π^2*T)); #(*fm^-3/GeV*) added a QC term for J/Psi
            n01+= QC*QC*degeneracy*(m^2* T /(2 *π^2)* ex* b2); #(*fm^-3/GeV*) QC term for J/Psi
            n11 +=QC*QC*degeneracy*(ex*m^2*(m*b1 + 2*T*b2 + m*b3)/(4*π^2*T));
               
        else
            ## loop over barions 
            b2 = besselkx(2,reducemass)
            b1 = besselk1x(reducemass)
            b3 = b1+4/(reducemass)*b2
            ex=exp(QC* μ - reducemass)
            density += degeneracy*(m^2* T /(2 *π^2)* ex* b2); #fm-3
            n10 += degeneracy*((ex*m^2*(m*b1 +2*T*b2 + m*b3))/(4*π^2*T)); #(*fm^-3/GeV*)
            n01+= QC*degeneracy*(m^2* T /(2 *π^2)* ex* b2); #(*fm^-3/GeV*)
            n11 +=QC*degeneracy*(ex*m^2*(m*b1 + 2*T*b2 + m*b3)/(4*π^2*T));
       
            
        end
    end
end
    return Thermodynamic(density*fmGeV3,(n10*fmGeV3,n01*fmGeV3),(n20*fmGeV3,n11*fmGeV3,n02*fmGeV3))
end

function hq_pressure(T,α;m=1.5)
    n = 1
    QC = 1
    density=zero(T)
    n10=zero(T)
    n01=zero(T)
    n20=zero(T)
    n02=zero(T)
    n11=zero(T)

    for i in eachindex(x) 
        ## loop over all particles 
        QC = round((x[i].Nc+x[i].Nac))
        m=x[i].Mass
        
        reducemass=m/T
        Spin=x[i].Spin
        degeneracy=x[i].Degeneracy
        
        if reducemass<500*one(T)
           
                b2 = besselkx(2,reducemass)
                b1 = besselkx(1,reducemass)
                
                b3 = b1+4/(reducemass)*b2
                ex=exp(QC* α - reducemass)
                #correction factor due to canonical ensemble 
                if QC == 1 
                    fact = besseli(1, ccbar/2)./besseli(0, ccbar/2)

                else 
                    fact = 1
                end   

                density += fact*QC*degeneracy*(m^2* T /(2 *π^2)* ex* b2); #fm-3
                n10 += fact*QC*degeneracy*((ex*m^2*(m*b1 +2*T*b2 + m*b3))/(4*π^2*T)); #(*fm^-3/GeV*)
                n01+= fact*QC*QC*degeneracy*(m^2* T /(2 *π^2)* ex* b2); #(*fm^-3/GeV*)
                n11 +=fact*QC*QC*degeneracy*(ex*m^2*(m*b1 + 2*T*b2 + m*b3)/(4*π^2*T));
            
        else 
                    density+= zero(T)
                    n01+= zero(T)
                    n10+= zero(T)
                    n11+= zero(T)

        end
    end

    return Thermodynamic(density*fmGeV3,(n10*fmGeV3,n01*fmGeV3),(n20*fmGeV3,n11*fmGeV3,n02*fmGeV3))
end




function invert_n_for_alpha(T, n; m=1.5, g=6)
    x = m/T
    A = g*m^2/(2π^2)

    K2 = besselkx(2, x)

    # For single-species Boltzmann with scaled Bessels:
    # n = A*T*exp(α-x)*K2  ->  exp(α-x) = n/(A*T*K2)
    # α = log(exp(α-x)) + x
    α = (n > 0) ? (log(n/(A*T*K2)) + x) : -Inf

    return α
end

#function hq_pressure(T, α; m=1.5, g=6)
#    x = m/T
#    A = g*m^2/(2π^2)
#
#    K0x = besselkx(0, x)
#    K1x = besselkx(1, x)
#    K2x = besselkx(2, x)
#    K3x = besselkx(3, x)
#    K4x = besselkx(4, x)
#
#    E = exp(α - x)            
#
#    P  = A*T^2*E*K2x
#    Pα = P
#    PT = A*E*( 2T*K2x + 0.5*m*(K1x + K3x) )  
#
#    Pαα = P
#    PTα = PT
#
#    PTT = A*E*( 2*K2x + x*(K1x + K3x) + (x^2/4)*(K0x + 2K2x + K4x) )
#
#    return Thermodynamic(P, (PT, Pα), (PTT, PTα, Pαα))
#end

