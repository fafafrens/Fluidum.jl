using Plots

data = hdf5_to_ObservableResult(pwd()*"/event_by_event_results_debug.h5")
glauber_vec = extract_glauber_multiplicity(data)
vn_vector = extract_vn(data)


spectra_pion1 = vn_vector[1,3,:,1,1]
spectra_pion2 = vn_vector[2,3,:,1,1]

spectra_D01 = vn_vector[1,3,:,1,2]
spectra_D02 = vn_vector[2,3,:,1,2]

sum(spectra_pion1)
sum(spectra_pion2)

sum(spectra_D01)
sum(spectra_D02)
plot(spectra_D01)
plot!(spectra_D02)

p = plot(pt_list, spectra_pion[1,:],yscale = :log10)
for i in 2:cc_ev_num
   plot!(pt_list, spectra_pion[i,:],yscale = :log10, legend=false)
end

spectra_D0 = vn_vector[sorted_indices,3,:,1,2]
plot!(pt_list, spectra_D0[1,:],yscale = :log10)
for i in 2:cc_ev_num
   plot!(pt_list, spectra_D0[i,:],yscale = :log10, legend=false)    
end
display(p)

# CENTRALITY CLASSES

"""
divides data into cc_fraction centrality classes based on final particle multiplicity
WIP ...
"""
function select_cc_events_multiplicity(data,cc_fraction)
    vn_vec = extract_vn(data)
    mult_vec = [total_M(result_single_event,species_list) for result_single_event in data]
    sorted_indices = sortperm(mult_vec,rev=true)
    Nev = length(mult_vec)
    cc_ev_num = div(cc_fraction*Nev,100)
    selected_data = [data[sorted_indices[1+i*cc_ev_num:cc_ev_num+i*cc_ev_num]] for i in 0:cc_fraction-1]
    return selected_data
end
data
data_chunks = select_cc_events_glauber(data,10)

q_vector_event_pt_dependent(data_chunks[1][1],species_list,[2,3])
q_vector_event_pt_dependent_im(data_chunks[1][1],species_list,[2,3])
multiplicity_event(data_chunks[1][1],species_list)
g = g_species_event_pt_dependent(data_chunks[1][1],species_list)
sum(g)
data_chunks[1]
q_vector_event_integrated(data_chunks[5][1],species_list,[2,3])

vns = [harmonic_coefficient(data_chunks[i],species_list,[2,3]) for i in 1:10]
vns[1]

using LaTeXStrings
default(lw = 2, size=(800,600),xtickfontsize=16,ytickfontsize=16,xlabelfontsize=16,ylabelfontsize=16,legendfontsize=16,grid=false,framestyle=:box)

plot(vns[1].vm_result[:,1,1], label = L"v_2\, \mathrm{\pi}", xlabel = L"p_T\, \mathrm{[GeV]}", ylabel = L"v_n")
plot!(vns[1].vm_result[:,1,2], label = L"v_3\, \mathrm{\pi}", xlabel = L"p_T\, \mathrm{[GeV]}", ylabel = L"v_n")
plot!(vns[1].vm_result[:,2,1], label = L"v_2\, \mathrm{D^0}", xlabel = L"p_T\, \mathrm{[GeV]}", ylabel = L"v_n")
plot!(vns[1].vm_result[:,2,2], label = L"v_3\, \mathrm{D^0}", xlabel = L"p_T\, \mathrm{[GeV]}", ylabel = L"v_n")
plot!(vns[1].vm_result_charged[:,1,1], label = L"v_2\, \mathrm{charged}", xlabel = L"p_T\, \mathrm{[GeV]}", ylabel = L"v_n")
plot!(vns[1].vm_result_charged[:,1,2], label = L"v_3\, \mathrm{charged}", xlabel = L"p_T\, \mathrm{[GeV]}", ylabel = L"v_n",legendtitle = L"\mathrm{0-10\% \,O-O}")


plot(vns[3].vm_result[:,1,1], label = L"v_2\, \mathrm{\pi}", xlabel = L"p_T\, \mathrm{[GeV]}", ylabel = L"v_n")
plot!(vns[3].vm_result[:,1,2], label = L"v_3\, \mathrm{\pi}", xlabel = L"p_T\, \mathrm{[GeV]}", ylabel = L"v_n")
plot!(vns[3].vm_result[:,2,1], label = L"v_2\, \mathrm{D^0}", xlabel = L"p_T\, \mathrm{[GeV]}", ylabel = L"v_n")
plot!(vns[3].vm_result[:,2,2], label = L"v_3\, \mathrm{D^0}", xlabel = L"p_T\, \mathrm{[GeV]}", ylabel = L"v_n")
plot!(vns[3].vm_result_charged[:,1,1], label = L"v_2\, \mathrm{charged}", xlabel = L"p_T\, \mathrm{[GeV]}", ylabel = L"v_n")
plot!(vns[3].vm_result_charged[:,1,2], label = L"v_3\, \mathrm{charged}", xlabel = L"p_T\, \mathrm{[GeV]}", ylabel = L"v_n",legendtitle = L"\mathrm{20-30\% \,O-O}")


spectra_event(data_chunks[1][1],species_list)
spectra(data_chunks[1],species_list)
plot(spectra(data_chunks[1],species_list)[:,1], xlabel = L"p_T\, \mathrm{[GeV]}", ylabel = L"\frac{dN}{dyp_Tdp_T}",yscale=:log10)
plot!(spectra(data_chunks[2],species_list)[:,1], xlabel = L"p_T\, \mathrm{[GeV]}", ylabel = L"\frac{dN}{dyp_Tdp_T}",yscale=:log10)
plot!(spectra(data_chunks[5],species_list)[:,1], xlabel = L"p_T\, \mathrm{[GeV]}", ylabel = L"\frac{dN}{dyp_Tdp_T}",yscale=:log10)

plot!(spectra(data_chunks[1],species_list)[:,2], xlabel = L"p_T\, \mathrm{[GeV]}", ylabel = L"\frac{dN}{dyp_Tdp_T}",yscale=:log10)
plot!(spectra(data_chunks[2],species_list)[:,2], xlabel = L"p_T\, \mathrm{[GeV]}", ylabel = L"\frac{dN}{dyp_Tdp_T}",yscale=:log10)
plot!(spectra(data_chunks[5],species_list)[:,2], xlabel = L"p_T\, \mathrm{[GeV]}", ylabel = L"\frac{dN}{dyp_Tdp_T}",yscale=:log10)

plot!(vns[1].vm_result[:,2,1])
plot!(vns[1].vm_result[:,2,2])

plot!(vns[1].vm_result_charged[:,1,1])
plot!(vns[1].vm_result_charged[:,1,2])

plot([vns[i].vm_result_integrated[1,1,1] for i in 1:10],xlabel = "Centrality class", label = L"v_2\, \mathrm{\pi}")
plot!([vns[i].vm_result_integrated[1,2,1] for i in 1:10],xlabel = "Centrality class", label = L"v_2\, \mathrm{D^0}")
plot!([vns[i].vm_result_charged_integrated[1,1,1] for i in 1:10],xlabel = "Centrality class", label = L"v_2\, \mathrm{charged}")

plot!([vns[i].vm_result_integrated[1,2,1] for i in 1:10])


