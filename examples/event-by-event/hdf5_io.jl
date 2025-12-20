using HDF5

struct ObservableResult
    glauber_multiplicity::Float64
    pt_list::Vector{Float64}
    vn::Array{Float64,4} # (cos, sin, denom) x length(pTlist) x length(wavenum_m) x length(species_list)
end

"""
    append_to_dataset!(file, dataset_name, new_data)

Append `new_data` to an existing HDF5 dataset, extending along the first dimension.
"""
function append_to_dataset!(file, dataset_name::String, new_data::AbstractArray)
    dset = file[dataset_name]
    curr_dims = size(dset)
    n_new = size(new_data, 1)
    
    new_dims = (curr_dims[1] + n_new, curr_dims[2:end]...)
    HDF5.set_extent_dims(dset, new_dims)
    
    slab_size = ntuple(length(curr_dims)) do I
        I == 1 ? (curr_dims[1] + 1:new_dims[1]) : (1:curr_dims[I])
    end
    
    dset[slab_size...] = new_data
end

"""
    create_extensible_dataset!(file, dataset_name, data)

Create a new HDF5 dataset that can be extended along the first dimension.
"""
function create_extensible_dataset!(file, dataset_name::String, data::AbstractArray)
    initial_size = size(data)
    max_size = (-1, initial_size[2:end]...)
    chunk = initial_size
    
    dspace = dataspace(initial_size; max_dims=max_size)
    dset = HDF5.create_dataset(file, dataset_name, Float64, dspace, chunk=chunk)
    
    indices = ntuple(i -> i == 1 ? (1:initial_size[1]) : Colon(), ndims(data))
    dset[indices...] = data
end

"""
    extract_glauber_multiplicity(data)

Extract glauber multiplicity values from ObservableResult vector.
"""
extract_glauber_multiplicity(data::Vector{ObservableResult}) = 
    [d.glauber_multiplicity for d in data]

"""
    extract_pt_list(data)

Extract pt_list values from ObservableResult vector.
"""
extract_pt_list(data::Vector{ObservableResult}) = 
    [d.pt_list[i] for d in data, i in 1:length(data[1].pt_list)]

"""
    extract_vn(data)

Extract vn values from ObservableResult vector.
"""
extract_vn(data::Vector{ObservableResult}) = 
    [d.vn[i,j,k,l] for d in data, i in 1:3, j in 1:size(data[1].vn,2), k in 1:size(data[1].vn,3), l in 1:size(data[1].vn,4)]

"""
    append_to_h5(filename, data)

Append ObservableResult data to an HDF5 file. Creates the file if it doesn't exist.
"""
function append_to_h5(filename, data::Vector{ObservableResult})
    glauber_data = extract_glauber_multiplicity(data)
    pt_data = extract_pt_list(data)
    vn_data = extract_vn(data)
    
    if isfile(filename)
        h5open(filename, "r+") do file
            append_to_dataset!(file, "glauber_multiplicity", glauber_data)
            append_to_dataset!(file, "pt_list", pt_data)
            append_to_dataset!(file, "vn", vn_data)
        end
    else
        h5open(filename, "w") do file
            create_extensible_dataset!(file, "glauber_multiplicity", glauber_data)
            create_extensible_dataset!(file, "pt_list", pt_data)
            create_extensible_dataset!(file, "vn", vn_data)
        end
    end
end
