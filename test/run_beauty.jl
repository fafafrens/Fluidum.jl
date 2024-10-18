using Fluidum
dic.D0
cd("C:\\Users\\feder\\.julia\\dev\\Fluidum\\test")
particle_list=pwd()*"/../src/kernels/particles_beauty_pythia8_list.txt"
const M=4.75
eos = Heavy_Quark(mass=M,hadron_list=HadronResonaceGas(name_file=particle_list,Maxmass=11.,Minmass=4.8))
@show eos.hadron_list
eos_hadrons=HadronResonaceGas(name_file=particle_list,Maxmass=11.,Minmass=4.9)
Fluidum.dic_NP
begin
cd("C:\\Users\\feder\\.julia\\dev\\Fluidum\\test")
particle_list=pwd()*"/../src/kernels/particles_beauty_pythia8_list.txt"
eos = Heavy_Quark(mass=M,hadron_list=HadronResonaceGas(name_file=particle_list,Maxmass=11.,Minmass=4.9))
#@show eos.hadron_list
eos_hadrons=HadronResonaceGas(name_file=particle_list,Maxmass=11.,Minmass=4.9)
#thermodynamic(0.6,3,eos_hadrons; ccbar = .7)
#eos.mass

#=results = Fluidum.runFluidum(eos,eos_HQ=eos_hadrons,
ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.0,M=M,
tau0=0.4,maxtime=10,
gridpoints=500,rmax=30,
temp_profile=:t_int, fug_profile=:fug, fug_pars=Fluidum.FugacityPars(rdrop=7,norm_coll=1.,dσ_QQdy=0.0296))
#[Fluidum.plot_field(results,i; tspan=(0.4,10),Δt=2,gridpoints=500,rmax=30,save=true) for i in Fluidum.names(Fluidum.HQ_viscous_1d())]
#Fluidum.plot_density(results,eos=eos_hadrons, tspan=(0.4,10),save=true)
#@show Fluidum.test_integral_cauchy(results,CartesianDiscretization(OriginInterval(500,30)).grid ,(0.4,10),eos;dt=0.5)
#[obs1=Fluidum.compute_observables(eos,Fluidum.dic_NP[particle],M=M,DsT=0.0,eos_HQ=eos_hadrons,Tfo=0.156,save = true,savepath="./HP_results/",fug_pars=Fluidum.FugacityPars(rdrop=7,norm_coll=1.,dσ_QQdy=0.0296),pt_min=0.01,rmax=30) for particle in keys(Fluidum.dic_NP)]
#[obs1=Fluidum.compute_observables(eos,Fluidum.dic_P_pythia[particle],M=M,DsT=0.0,eos_HQ=eos_hadrons,Tfo=0.156,save = true,savepath="./HP_results/",fug_pars=Fluidum.FugacityPars(rdrop=7,norm_coll=1.,dσ_QQdy=0.0296),pt_min=0.01,rmax=30) for particle in keys(Fluidum.dic_P_pythia)]
=#
#obs1 = Vector{Any}(undef,length(keys(Fluidum.dic_P_pythia)))
#=
obs2 = Vector{Any}(undef,length(keys(Fluidum.dic_P_pythia)))
i = 1
   for particle in keys(Fluidum.dic_P_pythia)
        obs2[i]=Fluidum.compute_observables(eos,Fluidum.dic_P_pythia[particle],M=M,DsT=0.0,eos_HQ=eos_hadrons,Tfo=0.156,save = true,savepath="./HP_results/low_fug/",fug_pars=Fluidum.FugacityPars(rdrop=7,norm_coll=1.,dσ_QQdy=0.0296),pt_min=0.01,rmax=30,ccbar=1.4)
        #obs1[i]=Fluidum.compute_observables(eos,Fluidum.dic_P_pythia[particle],M=M,DsT=0.0,eos_HQ=eos_hadrons,Tfo=0.156,save = true,savepath="./HP_results/high_fug/",fug_pars=Fluidum.FugacityPars(rdrop=7,norm_coll=1.,dσ_QQdy=0.0296),pt_min=0.01,rmax=30,ccbar=1.4)
        #@show obs1[i]
        #Fluidum.plot_spectra(obs1[i],obs2[i],norm_spectra=1, thermal = false, total = true, save = true,path="./plots/fill/")
        #mult_prompt_low+=obs1[i].yield_tot
        #mult_prompt_up+=obs2[i].yield_tot
        i+=1
    end
  
=#
#obs1_np = Vector{Any}(undef,length(keys(Fluidum.dic_NP)))
#obs2_np = Vector{Any}(undef,length(keys(Fluidum.dic_NP)))
i=1
for particle in keys(Fluidum.dic_NP)
    #obs1_np[i]=Fluidum.compute_observables(eos,Fluidum.dic_NP[particle],M=M,DsT=0.0,eos_HQ=eos_hadrons,Tfo=0.156,save = true,savepath="./HP_results/low_fug",fug_pars=Fluidum.FugacityPars(rdrop=7,norm_coll=1.,dσ_QQdy=0.0296),pt_min=0.01,rmax=30,ccbar=1.4)
    #obs2_np[i]=Fluidum.compute_observables(eos,Fluidum.dic_NP[particle],M=M,DsT=0.0,eos_HQ=eos_hadrons,Tfo=0.156,save = true,savepath="./HP_results/high_fug",fug_pars=Fluidum.FugacityPars(rdrop=7,norm_coll=1.,dσ_QQdy=0.0296),pt_min=0.01,rmax=30,ccbar=1.4)
    Fluidum.plot_spectra(obs1_np[i],obs2_np[i],norm_spectra1=mult_prompt_high/mult_nprompt_high,norm_spectra2=mult_prompt_low/mult_nprompt_low, thermal = false, total = true, save = true,path="./plots_HP/")
    i+=1
   
end

    

#obs2=Fluidum.compute_observables(eos,Fluidum.dic_P_pythia[:Upsilon],M=M,DsT=0.0,eos_HQ=eos_hadrons,Tfo=0.156,save = true,fug_pars=Fluidum.FugacityPars(rdrop=7,norm_coll=1.,dσ_QQdy=0.0296),pt_min=0.01,rmax=50)
#@show obs1.yield_th
#@show obs1.yield_tot

end
Fluidum.plot_spectra(obs1_np[2],obs2_np[2],norm_spectra1=mult_prompt_high/mult_nprompt_high,norm_spectra2=mult_prompt_low/mult_nprompt_low, thermal = false, total = true, save = true,path="./plots_HP/")

mult_prompt_low = 0.
mult_prompt_high = 0.
for i in eachindex(obs1)
    mult_prompt_low+=obs1[i].yield_tot 
    mult_prompt_high+=obs2[1].yield_tot 
end    
mult_prompt_high
mult_prompt_low
mult_nprompt_low = 0.
mult_nprompt_high = 0.

for i in eachindex(obs1_np)
    mult_nprompt_low+=obs1_np[i].yield_tot 
    mult_nprompt_high+=obs2_np[1].yield_tot 
end  
mult_nprompt_low 
mult_nprompt_high

length(obs1)
length(obs2)
length(Fluidum.dic_P_pythia)
i=1
for particle in keys(Fluidum.dic_P_pythia)
    Fluidum.plot_spectra(obs1[i],obs2[i], thermal = true, total = true, save = true,path="./plots_HP/")
    i+=1
end
#1.05e9 -> 20.77
#0.86e9 -> 20.57
obs1=Fluidum.compute_observables(eos,Fluidum.dic_NP[:D0],M=M,DsT=0.0,eos_HQ=eos_hadrons,Tfo=0.156,save = false,savepath="./HP_results/",fug_pars=Fluidum.FugacityPars(rdrop=7,norm_coll=1.,dσ_QQdy=0.0296),pt_min=0.01,rmax=30,ccbar=1.4)
Fluidum.plot_spectra(obs1,norm_spectra=mult_prompt/mult_NP, thermal = false, total = true, save = true)
[obs1[1],obs1[2],obs1[3],obs1[4]]
Fluidum.plot_int_yields([obs1[1],obs1[2],obs1[3],obs1[5]],[obs2[1],obs2[2],obs2[3],obs2[5]])

mult_prompt_low
mult_prompt_up
mult_prompt =0.63
mult_NP = 0.3346
Fluidum.plot_density(results,eos=eos_hadrons, tspan=(0.4,10),save=true)
[Fluidum.plot_field(results,i; tspan=(0.4,10),Δt=2,gridpoints=500,rmax=30,save=true) for i in Fluidum.names(Fluidum.HQ_viscous_1d())]
temp = [0 0; 0 0]
Fluidum.jgemvavx(temp,)
resultsFO = Fluidum.runFluidum_fo(eos,eos_HQ=eos_hadrons,
ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.0,M=M,
tau0=0.4,Tfo=0.156,
gridpoints=500,rmax=30,
temp_profile=:t_int, fug_profile=:fug, fug_pars=Fluidum.FugacityPars(dσ_QQdy=0.0296))
keys(Fluidum.dic_P_pythia)
resultsFO.fo.fields
x,phi=resultsFO.fo
phi.a[100][6]
plot!([phi.a[i][6] for i in 1:100])
obs=Fluidum.compute_observables(eos,0.139,eos_HQ=eos_hadrons,Tfo=0.156,save = true,fug_pars=Fluidum.FugacityPars(norm_coll=.5,dσ_QQdy=0.0296))
obs.yield_th
obs=Fluidum.compute_observables(eos,5.2798,eos_HQ=eos_hadrons,Tfo=0.156,save = true,fug_pars=Fluidum.FugacityPars(rdrop = 4,norm_coll=1,dσ_QQdy=0.0296))
obs = [Fluidum.compute_observables(eos,Fluidum.dic_P[key],M=M,DsT=0.0,eos_HQ=eos_hadrons,Tfo=0.156,save = true,fug_pars=Fluidum.FugacityPars(norm_coll=.5,dσ_QQdy=0.0296),pt_min=0.01) for key in keys(dic_P)]
obs.yield_th
obs.yield_tot

obs1.yield_tot
obs2.yield_tot

Fluidum.test_integral_cauchy(results,CartesianDiscretization(OriginInterval(500,30)).grid,(1,10),eos;dt=2.)
obs1=Fluidum.compute_observables(eos,Fluidum.dic_P_pythia[:Upsilon_2S_],M=M,DsT=0.0,eos_HQ=eos_hadrons,Tfo=0.156,save = true,fug_pars=Fluidum.FugacityPars(rdrop=6,norm_coll=.5,dσ_QQdy=0.0296),pt_min=0.01)
obs2=Fluidum.compute_observables(eos,Fluidum.dic_P_pythia[:Upsilon],M=M,DsT=0.0,eos_HQ=eos_hadrons,Tfo=0.156,save = true,fug_pars=Fluidum.FugacityPars(rdrop=6,norm_coll=.5,dσ_QQdy=0.0296),pt_min=0.01)


keys(dic_P)
dic.Ups10000zer
Fluidum.plot_spectra(obs,save=true,thermal=false,total=true)
    

hey=Fluidum.Observables(resultsFO.fo,Fluidum.dic.Ups10000zer,resultsFO.fluidproperties,0.1565,pt_min=0.,pt_max=10.0,step=100)
Fluidum.test_integral_cauchy(results,CartesianDiscretization(OriginInterval(500,30)).grid ,(0.4,10),eos;dt=0.5)
hey.spectra_th
Fluidum.names(Fluidum.HQ_viscous_1d())

keys(dic_P)

#charm routine
begin

cd("C:\\Users\\feder\\.julia\\dev\\Fluidum\\test")
particle_list=pwd()*"/../src/kernels/particles_beauty_pythia8_list.txt"
#const M=4.75
eos = Heavy_Quark()
eos_hadrons=HadronResonaceGas()

results = Fluidum.runFluidum(eos,eos_HQ=eos_hadrons,
ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.0,
tau0=0.4,maxtime=10,
gridpoints=500,rmax=30,
temp_profile=:t_int, fug_profile=:fug, fug_pars=Fluidum.FugacityPars())
#[Fluidum.plot_field(results,i; tspan=(0.4,10),Δt=2,gridpoints=500,rmax=30,save=true) for i in Fluidum.names(Fluidum.HQ_viscous_1d())]
#Fluidum.plot_density(results, tspan=(0.4,10),save=true)
Fluidum.test_integral_cauchy(results,CartesianDiscretization(OriginInterval(500,30)).grid ,(0.4,10),eos;dt=0.5)
resultsFO = Fluidum.runFluidum_fo(eos,eos_HQ=eos_hadrons,
ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.0,M=M,
tau0=0.4,Tfo=0.156,
gridpoints=500,rmax=30,
temp_profile=:t_int, fug_profile=:fug, fug_pars=Fluidum.FugacityPars())
resultsFO.fo.fields
x,phi=resultsFO.fo
phi.a[100][6]
plot([phi.a[i][6] for i in 1:100])
#=
obs=Fluidum.compute_observables(eos,0.139,eos_HQ=eos_hadrons,Tfo=0.156,save = true,fug_pars=Fluidum.FugacityPars(norm_coll=.5,dσ_QQdy=0.0296))
obs.yield_th
obs=Fluidum.compute_observables(eos,M,eos_HQ=eos_hadrons,Tfo=0.156,save = true,fug_pars=Fluidum.FugacityPars(norm_coll=.5,dσ_QQdy=0.0296))
obs = [Fluidum.compute_observables(eos,Fluidum.dic_P[key],M=M,DsT=0.0,eos_HQ=eos_hadrons,Tfo=0.156,save = true,fug_pars=Fluidum.FugacityPars(norm_coll=.5,dσ_QQdy=0.0296),pt_min=0.01) for key in keys(dic_P)]
obs1.yield_th
obs2.yield_th

obs1.yield_tot
obs2.yield_tot

Fluidum.test_integral_cauchy(results,CartesianDiscretization(OriginInterval(500,30)).grid,(1,10),eos;dt=2.)
obs1=Fluidum.compute_observables(eos,Fluidum.dic_P_pythia[:Upsilon_2S_],M=M,DsT=0.0,eos_HQ=eos_hadrons,Tfo=0.156,save = true,fug_pars=Fluidum.FugacityPars(rdrop=6,norm_coll=.5,dσ_QQdy=0.0296),pt_min=0.01)
obs2=Fluidum.compute_observables(eos,Fluidum.dic_P_pythia[:Upsilon],M=M,DsT=0.0,eos_HQ=eos_hadrons,Tfo=0.156,save = true,fug_pars=Fluidum.FugacityPars(rdrop=6,norm_coll=.5,dσ_QQdy=0.0296),pt_min=0.01)


keys(dic_P)
dic.Ups10000zer
Fluidum.plot_spectra(obs,save=true,thermal=false,total=true)
    

hey=Fluidum.Observables(resultsFO.fo,Fluidum.dic.Ups10000zer,resultsFO.fluidproperties,0.1565,pt_min=0.,pt_max=10.0,step=100)
Fluidum.test_integral_cauchy(results,CartesianDiscretization(OriginInterval(500,30)).grid ,(0.4,10),eos;dt=0.5)
hey.spectra_th
Fluidum.names(Fluidum.HQ_viscous_1d())

keys(dic_P)
=#
end

begin

#Fluidum.test_integral_cauchy(results,CartesianDiscretization(OriginInterval(500,30)).grid,(1,10),eos;dt=2.)
#obs1=Fluidum.compute_observables(eos,Fluidum.dic_charm[:Dc1865zer],DsT=0.0,eos_HQ=eos_hadrons,Tfo=0.156,save = true,fug_pars=Fluidum.FugacityPars(),pt_min=0.02, pt_max = 6.)
obs2=Fluidum.compute_observables(eos,Fluidum.dic_charm[:jp3096zer],DsT=0.2,eos_HQ=eos_hadrons,Tfo=0.156,save = true,fug_pars=Fluidum.FugacityPars(),pt_min=0.02, pt_max = 6.)
@show obs1.yield_th
@show obs2.yield_th

@show obs1.yield_tot
@show obs2.yield_tot
Fluidum.plot_spectra(obs2,save=true,thermal=false,total=true)
end
#Fluidum.plot_spectra(obs2,save=true,thermal=false,total=true)

#keys(Fluidum.dic_charm)
#Fluidum.dic_charm[:Dc1865zer].name

