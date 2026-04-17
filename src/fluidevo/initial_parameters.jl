const dσ_QQdy_dict=Dict(
    ("ALICE", 5.02, "NNPDF") => [0.03952, 0.03639, 0.04264],
    ("ALICE", 5.02, "CTEQ") => [0.0463, 0.0372, 0.05540],
    ("ALICE", 5.02, "exp_shadowing") => [0.0757, 0.0600, 0.0924],
    ("RHIC", 0.200, "CTEQ") => [0.005968, 0.0023561, 0.012979],
)

const σ_in_dict=Dict(
    5.02=>7.,
    0.200 =>4.23,
)

const radius_dict=Dict(
    "PbPb"=>6.62,
    "AuAu"=>7.,
)


Base.@kwdef struct Detector{S, T}
    det_name::S #detector name
    energy::T #colliding energy
    nuclei::S #colliding nuclei
    PDF::S = "CTEQ" #PDF used for the charm cross section
end

dσ_QQdy(det::Detector) = dσ_QQdy_dict[(det.det_name, det.energy, det.PDF)][1]
σ_in(det::Detector) = σ_in_dict[det.energy]
nucl_radius(det::Detector) = radius_dict[det.nuclei]


Base.@kwdef struct GridParameters{T}
    rmax::T = 20.
    gridpoints::Int = 300
end



Base.@kwdef struct ExpTail{T}
    is_tail::Bool = true
    offset_x::T = 0.01
    x_max_x::T = 8.
    offset_y::T = 0.
    x_max_y::T = 5.
    rdrop::T = 8.
end

Base.@kwdef struct Centrality
    cent1::Int
    cent2::Int
end



Base.@kwdef struct RunConfig{S, T}
    det::Detector{S,T} = Detector() 
    grid::GridParameters{T} = GridParameters()
    temp_norm::T 
    tspan::Tuple{T,T}
    tail_pars::ExpTail{T} = ExpTail()
    cent::Centrality
    Tfo::T
    particle_list::String = HRGc_list
end



