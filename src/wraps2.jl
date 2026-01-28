#initial conditions
#obj: you want temp(x) and fug(x) that can depend on different params

abstract type AbstractInitialCondition end

#structures
struct Gaussian_Intial_Condition{T,S} <:AbstractInitialCondition
    σ::T
    mean::S
end 

struct Trento_Intial_Condition{T} <:AbstractInitialCondition
    norm::T
end 

struct Step_Intial_Condition{T,S} <:AbstractInitialCondition
    norm::T
    r0::S
end 

struct pQCD_Initial_Condition{T,S,U} <: AbstractInitialCondition
    norm_coll::T
    σ_in::S
    dσ_QQdy::U
end

struct Observables{S,T,U,V,M,K,A,B,C,D}
    particle::particle_attribute{S,T,U,V}
    yield_th::Float64
    yield_tot::Float64
    pt_bins::M
    spectra_th::K
    spectra_tot::K
    fluid_properties::FluidProperties{A,B,C,D}
    Tfo::Float64
end

struct Trento_Entropy_Intial_Condition{T,S,V} <:AbstractInitialCondition
    norm::T
    cent1::S
    cent2::S
    file::V
end



#defaults

charm_pQCD(;norm_coll=1, σ_in=70,dσ_QQdy=0.463) = pQCD_Initial_Condition(norm_coll, σ_in, dσ_QQdy)
beauty_pQCD(;norm_coll=1, σ_in=70,dσ_QQdy=0.0296) = pQCD_Initial_Condition(norm_coll, σ_in, dσ_QQdy)

charm_viscous_diff_fluidproperties(;eos=FluiduMEoS(),ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.2,M=1.5)=Fluidum.FluidProperties(eos,QGPViscosity(ηs,Cs),SimpleBulkViscosity(ζs,Cζ),HQdiffusion(DsT,M))
beauty_viscous_diff_fluidproperties(;eos=FluiduMEoS(),ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.2,M=4.8)=Fluidum.FluidProperties(eos,QGPViscosity(ηs,Cs),SimpleBulkViscosity(ζs,Cζ),HQdiffusion(DsT,M))
charm_viscous_fluidproperties(;eos=FluiduMEoS(),ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.,M=1.5)=Fluidum.FluidProperties(eos,QGPViscosity(ηs,Cs),SimpleBulkViscosity(ζs,Cζ),ZeroDiffusion())
beauty_viscous_fluidproperties(;eos=FluiduMEoS(),ηs=0.2,Cs=0.2,ζs=0.02,Cζ=15.0,DsT=0.0,M=4.8)=Fluidum.FluidProperties(eos,QGPViscosity(ηs,Cs),SimpleBulkViscosity(ζs,Cζ),ZeroDiffusion())

#methods

function temperature(position,ini::T) where {T<:AbstractInitialCondition}
    throw("temperature not defined for the type $T")
end


function temperature(position,ini::T) where {T<:Gaussian_Intial_Condition}
    exp(-(position-ini.mean)^2/2/ini.σ)
end

function temperature(position, ini::T) where {T<:Trento_Entropy_Intial_Condition}
    ini.norm*Profiles(TabulatedTrento(ini.file), ini.cent1, ini.cent2,xmax=15)(position)
end

function temperature(position,ini::T) where {T<:Trento_Intial_Condition}
    file = SMatrix{128,2}([
    50.00000000000001      0.0018209316033560715;
    40.13143650955929      0.002287314501830327;
    36.22370269718774      0.0025666672190076565;
    33.85197508381405      0.0027801925144112284;
    32.17442046843765      0.002958501213041191;
    30.88695317916113      0.0031146035323238587;
    29.846813228576732      0.0032553981525620473;
    28.976101912244275      0.0033850473415507534;
    28.227889367808377      0.0035062888723503727;
    27.571836640690563      0.0036210400188860952;
    26.987263738170615      0.0037307097189299993;
    26.45948039809319      0.0038363741839850017;
    25.977699880177557      0.00393888236628606;
    25.533783626042695      0.004038922669536505;
    25.1214502248107      0.004137067001626551;
    24.735757492541484      0.004233800925137927;
    24.37275222405351      0.004329544943039982;
    24.029226685409654      0.004424669948781903;
    23.702545208452385      0.004519508710195961;
    23.390518088457675      0.0046143646140958386;
    23.09130816798676      0.004709518461464867;
    22.803360486342374      0.004805233871309577;
    22.5253485133512      0.004901761654669856;
    22.25613250936185      0.004999343445790516;
    21.9947268871234      0.005098214774445739;
    21.7402743486504      0.005198607727243464;
    21.49202518530012      0.005300753310087607;
    21.249320558052624      0.005404883589330761;
    21.01157887848169      0.005511233686552283;
    20.778284628782345      0.00562004367473069;
    20.54897911767795      0.005731560425767166;
    20.323252785662174      0.005846039445770992;
    20.100738759838087      0.005963746735837068;
    19.881107423896115      0.006084960709239124;
    19.664061818338272      0.006209974198835031;
    19.449333724030403      0.0063390965844902615;
    19.236680311503836      0.006472656075698634;
    19.025881261276133      0.006611002181426785;
    18.81673627838642      0.006754508408567798;
    18.609062938503044      0.006903575225698277;
    18.402694814224773      0.007058633345582235;
    18.19747983921038      0.0072201473714962406;
    17.99327887502761      0.0073886198744807525;
    17.789964451486647      0.007564595967655678;
    17.5874196560039      0.007748668459402197;
    17.385537151452116      0.00794148368010796;
    17.184218305167143      0.00814374809206292;
    16.98337241443518      0.008356235814353167;
    16.782916015985343      0.008579797213226999;
    16.58277226884503      0.008815368745925268;
    16.382870401448123      0.009063984267760608;
    16.183145215172043      0.009326788068594643;
    15.983536637563331      0.009605049946396013;
    15.783989319427194      0.0099001826934452;
    15.584452270733344      0.010213762449137643;
    15.384878530951637      0.01054755246666431;
    15.185224869995189      0.010903530967581763;
    14.985451516431812      0.011283923897289173;
    14.78552191003955      0.011691243593192626;
    14.585402476139066      0.012128334604166884;
    14.385062419444168      0.012598428201103282;
    14.184473535438558      0.013105207496589054;
    13.983610037518478      0.013652885572266972;
    13.782448398342408      0.014246299634576068;
    13.580967204004489      0.014891025017225655;
    13.379147019801893      0.01559351388677639;
    13.176970266500671      0.016361264860067407;
    12.974421106122437      0.017203031511764107;
    12.771485336378177      0.018129080077844733;
    12.568150292966473      0.019151509714544327;
    12.364404759034327      0.02028465268354947;
    12.16023888117001      0.021545577032058075;
    11.955644091360453      0.022954720930776586;
    11.750613034401866      0.024536695771354613;
    11.545139500302101      0.026321303576161234;
    11.339218361257572      0.02834482027433617;
    11.132845512827297      0.03065159153327519;
    10.926017818961816      0.03329595038695362;
    10.71873306057647      0.036344343013392705;
    10.510989887386897      0.039877220322074476;
    10.302787772749936      0.04398947595562329;
    10.094126971276062      0.04878663280546845;
    9.88500847899998      0.05290741553491948;
    9.675433995914558      0.06382043096916276;
    9.465405890689942      0.07566886675436467;
    9.254927167414849      0.08809407920801682;
    9.044001434210571      0.10041764298628211;
    8.83263287358063      0.1125333628751034;
    8.620826214370265      0.12472743316258399;
    8.408586705219927      0.13715197827720585;
    8.195920089406336      0.14981626707589135;
    7.982832580972949      0.16312695285820142;
    7.769330842059314      0.17807516460363312;
    7.555421961345767      0.1958061600879209;
    7.341113433536246      0.21659286470531688;
    7.126413139807729      0.23968795753695873;
    6.91132932916019      0.2643610315999489;
    6.695870600605634      0.29006615598420393;
    6.480045886139358      0.31625693553098944;
    6.263864434440456      0.3423661894267891;
    6.047335795252401      0.3678220592165537;
    5.830469804397872      0.39211192716655036;
    5.613276569385111      0.4148499045260683;
    5.395766455565949      0.43578934041641343;
    5.177950072808224      0.4548340537678695;
    4.959838262647786      0.472025449747684;
    4.741442085887394      0.48752728419112046;
    4.5227728106119045      0.501564161884873;
    4.303841900590962      0.5142427650724484;
    4.08466100404216      0.5255678062316925;
    3.865241942729129      0.5355676060932408;
    3.6455967013705193      0.5443159632357876;
    3.4257374173371216      0.551946554806182;
    3.2056763706155675      0.5586403933798669;
    2.985425974018184      0.5645148676176227;
    2.7649987636195617      0.5695693077981921;
    2.544407389401313      0.5737859644817784;
    2.3236646060873465      0.5772823821027658;
    2.1027832641527273      0.5803303062746771;
    1.8817763009898762      0.5831540678193577;
    1.660656732216491      0.5857555845753688;
    1.4394376431101061      0.5879948784171799;
    1.2181321801546996      0.5898227597330361;
    0.9967535426851927      0.5913414175938958;
    0.7753149746160524      0.5926022504085607;
    0.5538297562405278      0.5935347126428415;
    0.3323111960873209      0.5940928248768641;
    0.11077262282170842      0.5943142252688297       
    ])
    
    ini.norm*linear_interpolation(reverse(file[:,1]),reverse(file[:,2]); extrapolation_bc=Flat())(position).+0.0001
end


function ncoll(position,ini::T) where {T<:AbstractInitialCondition}
    throw("ncoll not defined for the type $T")
end

function ncoll(position,ini::T) where {T<:Step_Intial_Condition}
    if position<=ini.r0
        return ini.norm
    else 
        return 0.0001 #fm-2 since n0 [fm-2]
    end
    #you can also put them uniform
end

function dNdxdy(position,ini1::T,ini2::U) where {T<:Step_Intial_Condition,U<:pQCD_Initial_Condition}
    2*ini2.dσ_QQdy/ini2.σ_in*ncoll(position,ini1)
end

function dNdxdy(position,ini1::T,ini2::U) where {T<:AbstractInitialCondition,U<:AbstractInitialCondition}
    throw("dNdxdy not defined for the type $T,$U")
end

function dNdy(ini1::T,ini2::U;rmin=0.,rmax=30) where {T<:Step_Intial_Condition,U<:pQCD_Initial_Condition}
    quadgk(position->2π*position*dNdxdy(position,ini1,ini2),rmin,rmax)
end

function nhard(position,tau0,ini1::T,ini2::U,ini3::V) where {T<:AbstractInitialCondition,U<:AbstractInitialCondition,V<:AbstractInitialCondition}
    throw("nhard not defined for the type $T,$U,$V")
end 

function nhard(position,tau0,ini1::T,ini2::U,ini3::V) where {T<:Step_Intial_Condition,U<:pQCD_Initial_Condition,V<:Trento_Intial_Condition}
    ini2.norm_coll/tau0*dNdxdy(position,ini1,ini2)*(temperature(position,ini3)/temperature(0,ini3))^4 +0.001
    #set dNdxdy(0,ini1,ini2) if you want uniform ncoll and radial dependence given only by entropy profile
end 

function fugacity(position,ini::T) where {T<:AbstractInitialCondition}
    throw("fugacity not defined for the type $T")
end

function fugacity(position,tau0,ini1::T,ini2::U,ini3::V,eos) where {T<:Step_Intial_Condition,U<:pQCD_Initial_Condition,V<:Trento_Intial_Condition}
    if position<=ini1.r0
        return log(nhard(position,tau0,ini1,ini2,ini3)/(thermodynamic(temperature(position,ini3),0.0,eos).pressure)).+ 0.0001
        else return log(nhard(ini1.r0,tau0,ini1,ini2,ini3)/(thermodynamic(temperature(ini1.r0,ini3),0.0,eos).pressure)).+ 0.0001
        end
end 

function set_initial_conditions(ini1::T,ini2::U,ini3::V,eos_HQ,tau0;gridpoints=500,rmax=30) where {T<:Step_Intial_Condition,U<:pQCD_Initial_Condition,V<:Trento_Intial_Condition}
    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    oned_visc_hydro = Fluidum.HQ_viscous_1d()
    disc_fields = DiscreteFields(oned_visc_hydro,disc,Float64) 
    phi=set_array((x)->temperature(x,ini3),:temperature,disc_fields); #temperature initialization
    set_array!(phi,(x)->fugacity(x,tau0,ini1,ini2,ini3,eos_HQ),:α,disc_fields); #fugacity initialization
    NQQ̄,err= quadgk(x->2*pi*x*tau0*thermodynamic(Fluidum.temperature(x,ini3),Fluidum.fugacity(x,tau0,ini1,ini2,ini3,Fluidum.HadronResonaceGas()),eos_HQ).pressure,0,rmax,rtol=0.00001)
    @show NQQ̄
    #@show phi
return DiscreteFields(disc,disc_fields,phi)
end

function runFluidum_fo(ini1::T,ini2::U,ini3::V,fluidproperties::FluidProperties{A,B,C,D},eos_HQ,tau0;
    Tfo=0.156,maxtime=30, gridpoints=500,rmax=30) where {A,B,C,D,T<:Step_Intial_Condition,U<:pQCD_Initial_Condition,V<:Trento_Intial_Condition}
    fields=set_initial_conditions(ini1,ini2,ini3,eos_HQ,tau0;gridpoints=gridpoints,rmax=rmax)
    tspan = (tau0,maxtime)
    if fields.initial_field[1,1]<Tfo
        throw("Tfo = "*string(Tfo)*" MeV is larger than max temperature in the inital profile T0 = "*string(fields.initial_field[1,1]))
    end
    return freeze_out_routine(fields.discrete_field,Fluidum.matrix1d_visc_HQ!,fluidproperties,fields.initial_field,tspan,Tfo=Tfo)
end

function runFluidum(ini1::T,ini2::U,ini3::V,fluidproperties::FluidProperties{A,B,C,D},eos_HQ,tau0;
    maxtime=30, gridpoints=500,rmax=30) where {A,B,C,D,T<:Step_Intial_Condition,U<:pQCD_Initial_Condition,V<:Trento_Intial_Condition}
    fields=set_initial_conditions(ini1,ini2,ini3,eos_HQ,tau0;gridpoints=gridpoints,rmax=rmax)
    tspan = (tau0,maxtime)
    return oneshoot(fields.discrete_field,Fluidum.matrix1d_visc_HQ!,fluidproperties,fields.initial_field,tspan)

end

function Observables(fo::FreezeOutResult{M,N},part::particle_attribute{S,T,U,V},fluidproperties::FluidProperties{A,B,C,D},Tfo;
    pt_min=0,pt_max=10.0,step=100) where {M,N,S,T,U,V,A,B,C,D}
    spectra_th,err=spectra(fo,part,pt_max=pt_max,pt_min=pt_min,step=step,decays=false)
    spectra_tot,err=spectra(fo,part,pt_max=pt_max,pt_min=pt_min,step=step,decays=true)
    yield_th, err=multiplicity(fo,part,decays=false)
    yield_tot, err=multiplicity(fo,part,decays=true)
    pt_bins = range(pt_min,pt_max,step) 
    return Observables(part,yield_th,yield_tot,pt_bins,spectra_th,spectra_tot,fluidproperties,Tfo)
end

function Observables(fo::FreezeOutResult{M,N},m::Float64,fluidproperties::FluidProperties{A,B,C,D},Tfo;
    pt_min=0,pt_max=10.0,step=100,deg=1) where {M,N,A,B,C,D}
    spectra_th=getindex.(spectra_internal(m,fo,pt_max=pt_max,pt_min=pt_min,step=step,deg=deg)[:],1)
    spectra_tot=spectra_th
    yield_th, err=multiplicity(m,fo)
    yield_tot = yield_th
    pt_bins = range(pt_min,pt_max,step) 
    part = particle_attribute(string(m),m,deg,nothing,nothing)
    return Observables(part,yield_th,yield_tot,pt_bins,spectra_th,spectra_tot,fluidproperties,Tfo)
end

function save_observables(obj::Observables{S,T,U,V,M,K,A,B,C,D};path="./results/",overwrite=false) where {S,T,U,V,M,K,A,B,C,D}
    if isdir(path)
    else mkdir(path)
    end

    filename =  get_filename(obj,path=path)
    @show filename
    if isfile(filename)
        println("file already exists\n")
        if overwrite
            println("overwriting...\n")
        else
            return nothing
        end
    else
    open(filename, "w") do io 
        write(io, "# int_yield_th: "*string(obj.yield_th)*", int_yield_tot: "*string(obj.yield_tot)*"\n")
        write(io, "# pt\t th\t tot\t \n")
        writedlm(io, [obj.pt_bins obj.spectra_th obj.spectra_tot])
        close(io)
    end
end

#overwrite flag
end

function get_filename(obj::Observables{S,T,U,V,M,K,A,B,C,D};path="./results/",toplot=false) where {S,T,U,V,M,K,A,B,C,D}
    
    part = obj.particle.name
    ηs = obj.fluid_properties.shear.ηs
    ζs = obj.fluid_properties.bulk.ζs
    DsT = obj.fluid_properties.diffusion.DsT
   
    Tfo = obj.Tfo

    if toplot == true
        return path*string(part)*"_Tfo_"*string(Tfo)*"_ηs_"*string(ηs)*"_ζs_"*string(ζs)*"_DsT_"*string(DsT)*"_observables.pdf"
    else
        return path*string(part)*"_Tfo_"*string(Tfo)*"_ηs_"*string(ηs)*"_ζs_"*string(ζs)*"_DsT_"*string(DsT)*"_observables.txt"
    end
end


