
function total_M(result_single_event,species_list)
    M = 0.
    pTlist = result_single_event.pTlist
    vn = result_single_event.vn
    for k in eachindex(species_list)
        for pt_idx in eachindex(pTlist)
            M+=vn[3,pt_idx,1,k]
        end
    end
    return M
    
end


function M_species(result_single_event,species_list)
    M = zeros(length(species_list))
    pTlist = result_single_event.pTlist
    vn = result_single_event.vn
    for k in eachindex(species_list)
        for pt_idx in eachindex(pTlist)
            M[k]+=vn[3,pt_idx,1,k]
        end
    end
    return M
    
end

function M_species_ptbin(result_single_event,species_list)
    pTlist = result_single_event.pTlist
    vn = result_single_event.vn
    M = zeros(length(pTlist),length(species_list))
    for k in eachindex(species_list)
        for pt_idx in eachindex(pTlist)
            M[pt_idx,k]+=vn[3,pt_idx,1,k]
        end
    end  
    return M
end


function g_species_ptbin(result_single_event,species_list)
    return M_species_ptbin(result_single_event,species_list)./total_M(result_single_event,species_list)
end