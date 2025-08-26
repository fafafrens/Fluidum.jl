


struct detector{S, T}
    name::S #name of the detector
    r::T #nuclear radius
    σ_in::T #pp inelastic cross section
    dσ_QQdy::T #heavy quark cross section from FONLL
    nucl::S #colliding nuclei
end

struct InitialParameters{S}
    norm::S
    tau0::S
    tau_fs::S
    rdrop::S
    InitialParameters(norm::S, tau0::S, tau_fs::S, rdrop::S) where S = new{S}(norm, tau0, tau_fs, rdrop)
    InitialParameters(norm::S, tau0::S, rdrop::S) where S = new{S}(norm, tau0, tau0, rdrop)
end


struct GridParameters{S}
    rmax::S
    gridpoints::S
end


# Name	  radius      σ_in	  dσ_QQdy       nucl		
#ALICE       6.62     7.00         0.0757	Pb_Pb
#RHIC        7        4.23         0.005968	Au_Au   
#ALICE1      6.62     7.00         0.0463	Pb_Pb



#=struct norm_struct{S, T}
    cent_class::S #centrality class
    pi_mult::T #pion multiplicity
    norm::T
end=#




