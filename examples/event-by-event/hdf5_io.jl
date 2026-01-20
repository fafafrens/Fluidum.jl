
struct ObservableResult
    glauber_multiplicity::Float64
    vn::Array{Float64,4} # (cos, sin, denom) x length(pTlist) x length(wavenum_m) x length(species_list)
end

null_observable() = ObservableResult(0.0, Array{Float64}(undef, 3, 0, 0, 0))
function null_observable(wavenum,species_list)
    pt_length_max = maximum(length.(pt_list.(species_list)))
    ObservableResult(0.0, zeros(3, pt_length_max, length(wavenum), length(species_list)))
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
    append_to_jld2(filename, data)

Append ObservableResult data to a JLD2 file. Creates the file if it doesn't exist.
"""
function append_to_jld2(filename, data::Vector{ObservableResult})
    glauber_data = extract_glauber_multiplicity(data)
    pt_data = extract_pt_list(data)
    vn_data = extract_vn(data)

    if isfile(filename)
        jldopen(filename, "a") do file
            file["glauber_multiplicity"] = vcat(file["glauber_multiplicity"], glauber_data)
            file["pt_list"] = vcat(file["pt_list"], pt_data)
            file["vn"] = vcat(file["vn"], vn_data)
        end
    else
        jldsave(filename; glauber_multiplicity=glauber_data, pt_list=pt_data, vn=vn_data)
    end
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
    extract_vn(data)

Extract vn values from ObservableResult vector.
"""
extract_vn(data::Vector{ObservableResult}) = 
    [d.vn[i,j,k,l] for d in data, i in 1:3, j in 1:size(data[1].vn,2), k in 1:size(data[1].vn,3), l in 1:size(data[1].vn,4)]


"""
    append_to_h5(filename, data, metadata)

Append ObservableResult data to an HDF5 file. Creates the file if it doesn't exist. If the file doesn't exist, also writes metadata as attributes.
"""
function append_to_h5(filename, data::Vector{ObservableResult}, metadata::NamedTuple)
    glauber_data = extract_glauber_multiplicity(data)
   # pt_data = extract_pt_list(data)
    vn_data = extract_vn(data)
    
    if isfile(filename)
        meta2 = read_metadata(filename)
        if check_metadata(metadata, meta2)
            h5open(filename, "r+") do file
            append_to_dataset!(file, "glauber_multiplicity", glauber_data)
            #append_to_dataset!(file, "pt_list", pt_data)
            append_to_dataset!(file, "vn", vn_data)
            end
        else
            new_filename = replace(filename, ".h5" => "_new.h5")
            h5open(new_filename, "w") do file
            append_metadata!(file, metadata)
            create_extensible_dataset!(file, "glauber_multiplicity", glauber_data)
            #create_extensible_dataset!(file, "pt_list", pt_data)
            create_extensible_dataset!(file, "vn", vn_data)
            end
        end
        
    else
        h5open(filename, "w") do file
            append_metadata!(file, metadata)
            create_extensible_dataset!(file, "glauber_multiplicity", glauber_data)
            #create_extensible_dataset!(file, "pt_list", pt_data)
            create_extensible_dataset!(file, "vn", vn_data)
        end
    end
end

function append_metadata!(file, metadata)
        g = create_group(file, "metadata")
        attrs = attributes(g)

        for (k, v) in pairs(metadata)
            attrs[string(k)] = v
        
    end
end

function read_metadata(filename)
    h5open(filename, "r") do file
        attrs = attributes(file["metadata"])
        meta = (; (Symbol(k) => read(attrs[k]) for k in keys(attrs))...)
        return meta
    end
end

"""
    check_metadata(meta1, meta2)
"""
function check_metadata(meta1, meta2)
    for (k, v) in pairs(meta1)
        if haskey(meta2, k)
            if meta2[k] != v
                @warn "Metadata mismatch for key $k: $v != $(meta2[k]). Creating new file..."
                return false
            end
        else
            @error "Metadata key $k not found in second metadata"
            return false
        end
    end
    return true
end




function hdf5_to_ObservableResult(h5_file::AbstractString)
    isfile(h5_file) || error("HDF5 file not found: $h5_file")

    glauber_data, vn_data = h5open(h5_file, "r") do file
        (read(file["glauber_multiplicity"]), read(file["vn"]))
    end
        
        Nev = size(glauber_data, 1)
        return [ObservableResult(glauber_data[i], vn_data[i, :, :, :, :]) for i in 1:Nev]
end



