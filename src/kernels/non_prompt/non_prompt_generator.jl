# Import necessary libraries
using DelimitedFiles

# Function to read the data from file
function read_data(filename)
    return readdlm(filename, '\t', Float64, skipstart=1)
end

# Function to write the resulting data to a file
using Printf

function write_output_file(output_file, header, data)
    open(output_file, "w") do f
        # Write the header (assuming it's a comma or tab-separated string)
        write(f, header * "\n")
        
        # Write the data with fixed-width columns (10 characters for each column)
        for row in range(1,size(data)[1])
           @printf(f, "%.5f\t%.5f\t", data[row,1], data[row,2])
           # @show row
            for i in range(3,10)
                # Print each item formatted as a string with width 10
                @printf(f, "\t %.8e", data[row,i])
                if i != length(row)
                    write(f, "\t ")  # Add a space between columns
                end
            end
            write(f, "\n")  # Newline after each row
        end
    end
end





# Function to compute the difference between the two files and output the result
function compute_difference(file1, file2, output_file)
    # Read data from both files
    data1 = read_data(file1)
    data2 = read_data(file2)

    # Ensure the files have the same number of rows and columns
    @assert size(data1) == size(data2) "Files must have the same structure and size"

    # Prepare the result array
    result = similar(data1)
    
    # Loop over the rows
    for i in eachindex(data1[:,1])
        # First two columns (pT and Ur) remain unchanged
        result[i, 1] = data1[i, 1]  # pT
        result[i, 2] = data1[i, 2]  # Ur

        # Columns 3 to 10 will be the difference
        for j in 3:10
            result[i, j] = data1[i, j] - data2[i, j]
        end
    end

    # Define the header to be written in the output file
    header = "#     1:pT [GeV]    2:Ur    3:Keq 1    4:Keq 2    5:Kshear 1    6:Kshear 2    7:Kshear 3    8:Kshear 4    9:Kbulk 1    10:Kbulk 2\n"
    
    # Write the output to a file
    write_output_file(output_file, header, result)
end

# File paths (replace with actual file paths)
T = 0.1565
fluidum_dir="C:/Users/feder/.julia/dev/Fluidum/"
for i in ["D_s+","D0","Jpsi"]
    cd(fluidum_dir*"src/kernels/non_prompt/")
    file1 = i*"_total_T"*string(T)*"_Kj.out"
    file2 = i*"_thermal_T"*string(T)*"_Kj.out"
    output_file = i*"_nonprompt_T"*string(T)*"_Kj.out"
    #@assert isfile(output_file) "Output file already exists"

# Call the function to compute the difference and write to output
    compute_difference(file1, file2, output_file)
end
