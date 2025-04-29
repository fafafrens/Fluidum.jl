


struct detector{S, T}
    name::S #name of the detector
    r::T #nuclear radius
    σ_in::T #pp inelastic cross section
    dσ_QQdy::T #heavy quark cross section from FONLL
    nucl::S #colliding nuclei
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




