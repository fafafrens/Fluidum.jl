struct TrenTo_IC{S,T,U} #struct for all physical trento parameters
    projectile::T #projectile string
    p::S #reduced thickness
    k::S #fluctuations
    w::S #nucleon width
    v::S #constituent width
    m::U #number of constituents (INT)
    d::S #nucleon min dist
    x::S #cross section
end

struct TrenTo_IC_stats{T,U} #struct for statistics of IC
    NumEvents::T #min bias events for bins
    CentralityRange::U #list of centralities
    NumEventsAve::T #number of events for getting the profiles
end


#function that returns initial profiles
function get_initial_profile(ic::TrenTo_IC{T,S,U};path="./initial_conditions/",statistics=TrenTo_IC_stats(1e6,1:1:100,1e6)) where {T,S,U}
    flagbins, flagBG = check_for_files_bg(ic;path=path)
    nameBG, namebins = assembleTrentoName(ic;path=path)

    
    #run(`wsl`)
    #check if bins are available
    if flagbins
        bins=readdlm(namebins)
    else
        bins=get_entro_bins(ic,statistics;path=path)
    end

    #check if profiles are available
    if flagBG
        initial_profiles=readdlm(nameBG)
    else
        initial_profiles= run_average_eff(ic,bins;path=path,statistics=statistics)
    end

    return initial_profiles

end

#run events, sort them in the bins and take the event&angle average over them
function run_average_eff(ic::TrenTo_IC{T,S,U},bins;path="./initial_conditions/",statistics=TrenTo_IC_stats(1e6,1:1:100,1e6)) where {T,S,U}
    rList=collect(0:0.1:10)
    rLen=length(rList)
    centraLen=length(statistics.CentralityRange)
    profileMat=zeros(Float64,centraLen+1,rLen)
    profileMat[1,:]=rList
    indexList=zeros(Int,centraLen)
    trento_grid=-9.9:0.2:9.9
    tLen=length(trento_grid)
    ave_temp_mat=zeros(centraLen,tLen,tLen)

    for i in 1:1:statistics.NumEventsAve
        event_index=run_one_event_average(ic,bins)[1]
        indexList[event_index]=indexList[event_index]+1
        fid=h5open("temp.hdf","r")
        temp=read(fid)["event_0"]
        CoM = SciPy.ndimage.measurements.center_of_mass(temp) #get center of mass of profile
        shiftedTemp =SciPy.ndimage.shift(temp,(CoM[1]-50,CoM[2]-50)) # shift to center, change for different grid!
        phiR=rand(0:360) #rotation angle
        rotatedTemp = SciPy.ndimage.rotate(shiftedTemp, phiR, reshape=false)
        ave_temp_mat[event_index,:,:]=ave_temp_mat[event_index,:,:]+rotatedTemp
        close(fid)
        #cd("C://Users//anki1//Documents//PostDoc//Duke//Code//trento-master//build//src//")
        #run(`rm temp.hdf`)
        rm("temp.hdf")
    end
     #add check for additional reruns
     @show indexList
      @show sum(indexList)
    for i in 1:1:centraLen
        tempProfAve= ave_temp_mat[i,:,:] ./ indexList[i]
        interpolated_profile=Spline2D(trento_grid,trento_grid,tempProfAve)
        angle_average(r)=quadgk(phi->interpolated_profile(r*cos(phi),r*sin(phi)),0,2*pi)
        evaluated_profile=getfield.(angle_average.(rList),1)./(2*pi)
        profileMat[i+1,:]=evaluated_profile
    end
    fileNameBG, fileNameBins = assembleTrentoName(ic;path=path)
    writedlm(fileNameBG,profileMat)
    return profileMat
end

#gets name for file with specified trento params
function assembleTrentoName(ic::TrenTo_IC{T,S,U};path="./initial_conditions/") where {T,S,U}
    projectile=ic.projectile
    p=string(Int(10*ic.p))
    k=string(Int(10*ic.k))
    w=string(Int(10*ic.w))
    v=string(Int(10*ic.v))
    m=string(ic.m)
    d=string(Int(10*ic.d))
    x=string(Int(100*ic.x))
    fileNameBG=path*"trento_bg_"*projectile*"_x_"*x*"_p_"*p*"_k_"*k*"_w_"*w*"_v_"*v*"_m_"*m*"_d_"*d*".txt"
    fileNamebins=path*"trento_bins_"*projectile*"_x_"*x*"_p_"*p*"_k_"*k*"_w_"*w*"_v_"*v*"_m_"*m*"_d_"*d*".txt"
    return (fileNameBG,fileNamebins)
end

#returns bool to check if IC with specified params are available
function check_for_files_bg(ic::TrenTo_IC{T,S,U};path="./initial_conditions/") where {T,S,U}
    (fileNameBG,fileNamebins)=assembleTrentoName(ic;path=path)
     return (isfile(fileNamebins),isfile(fileNameBG))
end

function get_trento_entropy(ic::TrenTo_IC{S,T,U}) where {S,T,U}
    base1="wsl" 
    base2="C:/Users/anki1/Documents/PostDoc/Duke/Code/trento-master/build/src"
    base3="&&" 
    base4="trento"
    projectileString=ic.projectile
    k="-k"
    kval=string(ic.k)
    p="-p"
    pval=string(ic.p)
    x="-x"
    xval=string(ic.x);
    #mycmd1=`$base1 $base2`
    #run(mycmd1)
    #@show pwd()
    cd("C://Users//anki1//Documents//PostDoc//Duke//Code//trento-master//build//src//")
    #@show pwd()
    #run(`C:/Program Files/Git/bin/bash.exe -c "cd //mnt//c//Users//anki1//Documents//PostDoc//Duke//Code//trento-master//build//src//"`)#cd("//mnt//c//Users//anki1//Documents//PostDoc//Duke//Code//trento-master//build//src//")
    mycommand=`$base1 $base4 $projectileString $projectileString $k $kval $x $xval $p $pval`
    #run(`wsl trento Pb Pb`)
    #mycommand=`$base $projectileString $projectileString $optionsString`
    #@show mycommand
    trentoOutVec=parse.(Float64,filter!(e->e !="",split(readchomp(mycommand)," ")))
    return trentoOutVec[4]
end

#get entropy bins that define the centralities
function get_entro_bins(ic::TrenTo_IC{T,S,U},stats::TrenTo_IC_stats{A,B};path="./initial_conditions/") where {T,S,U,A,B}
    NumEvs=stats.NumEvents
    NumEvsPerBin=round(Int, NumEvs / length(stats.CentralityRange))
    entroList=Vector{Float64}(undef,NumEvs) 
    for i in 1:NumEvs
        tempEntro=get_trento_entropy(ic)
        entroList[i]=tempEntro
    end
    sort!(entroList,rev=true)
    bins=entroList[begin:NumEvsPerBin:end]
    nameBG, namebins = assembleTrentoName(ic;path=path)
    writedlm(namebins,bins)
    return bins
end

#generate one event and safe the profile & get the centrality index
function run_one_event_average(ic::TrenTo_IC{S,T,U},bins) where {S,T,U}
    base1="wsl" 
    base2="C://Users//anki1//Documents//PostDoc//Duke//Code//trento-master//build//src"
    base3="&&" 
    base4="trento"
    projectileString=ic.projectile
    k="-k"
    kval=string(ic.k)
    p="-p"
    pval=string(ic.p)
    x="-x"
    xval=string(ic.x)
    outbase="-o"
    outval="temp.hdf"
    #mycmd1=`$base1 $base2`
    #run(mycmd1)
    cd("C://Users//anki1//Documents//PostDoc//Duke//Code//trento-master//build//src//")
    mycommand=`$base1 $base4 $projectileString $projectileString $k $kval $x $xval $p $pval $outbase $outval`
    trentoOutVec=parse.(Float64,filter!(e->e !="",split(readchomp(mycommand)," ")))
    return get_bin_index(trentoOutVec[4],bins)
end

#return the index corresponding to the entropy/centrality bin of the event
function get_bin_index(entropy,bins)
    ind=findfirst(t -> t < entropy, bins)
    if isnothing(ind)
        return length(bins)
    else
        if ind[1]==1
            return ind[1]
        else
            return ind[1]-1
        end
    end
end