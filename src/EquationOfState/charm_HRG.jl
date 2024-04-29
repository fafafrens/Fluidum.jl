

struct HadronResonaceGas{T} <:EquationOfState
    particle_list::T
end

# Define a particle list, such as an array of particles (type T): particle_list = ["proton", "neutron", "pion"]
# Create an instance of HadronResonanceGas with the particle list: hadron_gas = HadronResonanceGas(particle_list)

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




#read in resonaces
function HadronResonaceGas(;name_file=string(@__DIR__,"/OpenCharmParticleList_corrJS.txt"),Maxmass=4,Minmass=1.0,condition=x->true)
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


#read in resonances
function HadronResonaceGasNew(;name_file="src/EquationofState/PDG2016Plus_massorder.dat",Maxmass=3.2,Minmass=1.8,condition=x->true)
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




function readresonancelist(;name_file=string(@__DIR__,"/OpenCharmParticleList_corrJS.txt"),Maxmass=4,Minmass=1.0,condition=x->true)
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

# #read in resonaces
# function HadronResonaceGas_ccbar(filterlist,ccbar;name_file=string(@__DIR__,"/OpenCharmParticleList_corrJS.txt"),Maxmass=4,Minmass=1.0,condition=x->true)
#     HadronResonaceGas_ccbar(filterlist,ccbar)

# end

#pretty printing
function Base.show(io::IO, z::HadronResonaceGas)
    min ,max = extrema(z.particle_list.Mass)
    print(io,"Hadron Resonace gas: ",length(z.particle_list)," particles with mass ⊆ ",min,"..",max," GeV" )
    for part in z.particle_list
        print(io,part.Name," mass=",part.Mass,"\t","spin=",part.Spin," QC=",(part.Nac+part.Nc)," deg=",(part.Degeneracy),"\n")
    end 
end


function Base.show(io::IO, ::MIME"text/plain", z::HadronResonaceGas) #where{T}
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


#show(HadronResonaceGas())
#condition to exclude particles, which are already included in Walecka model (also include neutrons?)
waleckacondition(x)=x.Name != "f00600zer"&&x.Name !="om0782zer"&& x.Name != "pr0938plu"  && x.Name !="pr0938plb" && x.Name !="ne0939zer" && x.Name != "ne0939zrb"#remove neutrons
#prova=HadronResonaceGas(name_file="EquationofState/particles.data",condition=waleckacondition)



@inline Base.getindex(elm::HadronResonaceGas, i::Int)= elm.particle_list[i] #to access the elements of particle_list
@inline Base.eachindex(elm::HadronResonaceGas)=Base.eachindex(elm.particle_list)

@inline Base.getindex(elm::HadronResonaceGasNew, i::Int)= elm.particle_list[i]
@inline Base.eachindex(elm::HadronResonaceGasNew)=Base.eachindex(elm.particle_list)

@inline Base.getindex(elm::HadronResonaceGas_ccbar, i::Int)= elm.particle_list[i]
@inline Base.eachindex(elm::HadronResonaceGas_ccbar)=Base.eachindex(elm.particle_list)


"""
    thermodynamic(T,μ,x::HadronResonaceGas{L})  where {L}

TBW test test hello
"""
function thermodynamic(T,μ,x::HadronResonaceGas{L}; ccbar = 2.76)  where {L}

    if isless(T,zero(T))
        return Thermodynamic(zero(T),(zero(T),zero(T)),(zero(T),zero(T),zero(T)))
    end

    density=zero(T)
    n10=zero(T)
    n01=zero(T)
    n20=zero(T)
    n02=zero(T)
    n11=zero(T)

    for i in eachindex(x) 
        ## loop over all particles 
        #QB = round((x[i].Nq+x[i].Ns+x[i].Nc-x[i].Naq-x[i].Nas-x[i].Nac)/3.0)
        #QS = round(x[i].Nas-x[i].Ns)
        QC = round((x[i].Nc+x[i].Nac))
        m=x[i].Mass
        
        reducemass=m/T
        Spin=x[i].Spin
        degeneracy=x[i].Degeneracy
        
        if reducemass<500*one(T)
           
                b2 = besselkx(2,reducemass)
                b1 = besselk1x(reducemass)
                
                
                #b2 = besselk(2,m/T)
                #b1 = besselk(1,m/T)
                b3 = b1+4/(reducemass)*b2
                ex=exp(QC* μ - reducemass)
                #equation 7, hq new
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
                    #p20+=zero(T)
                    n11+= zero(T)
                    #p02+= zero(T)

        end
    end

    return Thermodynamic(density*fmGeV3,(n10*fmGeV3,n01*fmGeV3),(n20*fmGeV3,n11*fmGeV3,n02*fmGeV3))
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

    density=zero(T)
    n10=zero(T)
    n01=zero(T)
    n20=zero(T)
    n02=zero(T)
    n11=zero(T)

    for i in eachindex(x) 
        ## loop over all particles 
        QB = x[i].Baryon
        #QS = round(x[i].Nas-x[i].Ns)
        QC = round((x[i].Charm))
        m=x[i].Mass
        reducemass=m/T
        degeneracy=x[i].Degeneracy*2
        
    if reducemass<1000*one(T)
        if QB==0
          
            ## loop over mesons
            b2 = besselkx(2,reducemass)
            b1 = besselk1x(reducemass)
            
            #b2 = besselk(2,m/T)
            #b1 = besselk(1,m/T)
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
            #b2 = besselk(2,m/T)
            #b1 = besselk(1,m/T)
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


"""
    thermodynamic(T,μ,x::HadronResonaceGas_ccbar{L})  where {L}

Add correction factor
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

    for i in eachindex(x) 
        ## loop over all particles 
        QC = round((x[i].Nc+x[i].Nac))
        m=x[i].Mass
        
        reducemass=m/T
        Spin=x[i].Spin
        degeneracy=x[i].Degeneracy
        
        if reducemass<500*one(T)
           
                b2 = besselkx(2,reducemass)
                b1 = besselk1x(reducemass)
                
                b3 = b1+4/(reducemass)*b2
                ex=exp(QC* μ - reducemass)
                #equation 7, hq new
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


function hq_pressure(T,μ;m=1.5)
    n = 1
    QC = 1
    density=zero(T)
    n10=zero(T)
    n01=zero(T)
    n20=zero(T)
    n02=zero(T)
    n11=zero(T)


    reducemass=m/T
    degeneracy= 6

    ## loop over barions 
    b2 = besselkx(2,reducemass)
    b1 = besselk1x(reducemass)
    #b2 = besselk(2,m/T)
    #b1 = besselk(1,m/T)
    b3 = b1+4/(reducemass)*b2
    ex=exp(QC* μ - reducemass)
    density += degeneracy*(m^2* T /(2 *π^2)* ex* b2); #fm-3
    n10 += degeneracy*((ex*m^2*(m*b1 +2*T*b2 + m*b3))/(4*π^2*T)); #(*fm^-3/GeV*)
    n01+= degeneracy*(m^2* T /(2 *π^2)* ex* b2); #(*fm^-3/GeV*)
    n11 +=degeneracy*(ex*m^2*(m*b1 + 2*T*b2 + m*b3)/(4*π^2*T));

   return Thermodynamic(density,(n10,n01),(n20,n11,n02))
end


