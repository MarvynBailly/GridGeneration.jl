"""
    readTurtleFields(filename::String)

Reads binary data from a "Turtle" format file, correctly handling row-major
data layout for all fields.

This function translates a Python script for reading structured grid data
from a custom binary format. It handles file endianness, reads a header,
and then iteratively reads different data fields (scalars, vectors, tensors)
defined on the grid blocks.

# Arguments
- `filename::String`: The path to the turtle file.

# Returns
A `NamedTuple` containing:
- `data::Dict{String, Any}`: A dictionary mapping variable names to their data arrays.
- `datatype::Dict{String, Int}`: A dictionary mapping variable names to their type identifiers.
- `gridfile::String`: The name of the associated grid file as stored in the header.

Returns `nothing` if a critical error occurs (e.g., wrong file format).
"""
function readTurtleFields(filename::String)
    # --- Parameters from Turtle_parameters.h ---
    MAGICNUMBER = 751123
    IO_NAME_LENGTH = 200

    # --- Parameters from VariableList.h ---
    CV_SCALAR = 1
    CV_VECTOR = 2
    CV_TENSOR = 3
    FA_SCALAR = 11
    FA_VECTOR = 12
    FA_TENSOR = 13
    NO_SCALAR = 21
    NO_VECTOR = 22
    NO_TENSOR = 23

    vartypeDict = Dict(
        "CV_SCALAR" => CV_SCALAR, "CV_VECTOR" => CV_VECTOR, "CV_TENSOR" => CV_TENSOR,
        "FA_SCALAR" => FA_SCALAR, "FA_VECTOR" => FA_VECTOR, "FA_TENSOR" => FA_TENSOR,
        "NO_SCALAR" => NO_SCALAR, "NO_VECTOR" => NO_VECTOR, "NO_TENSOR" => NO_TENSOR,
    )
    typekeyDict = Dict(v => k for (k, v) in vartypeDict)

    println("Reading file $filename")

    return open(filename, "r") do fid
        # --- Endianness Check ---
        magicno = read(fid, Int32)
        swap_endian = false
        if magicno != MAGICNUMBER
            magicno = bswap(magicno)
            if magicno == MAGICNUMBER
                swap_endian = true
                println("File is non-native endian. Will swap byte order.")
            else
                @error "Cannot find correct endian. Terminating reading of grid..."
                return nothing
            end
        end

        # Helper function to read a single value and swap if needed
        read_scalar(io, type) = swap_endian ? bswap(read(io, type)) : read(io, type)

        # Helper function to read an array and swap all elements if needed
        function read_array!(io, arr::AbstractArray)
            read!(io, arr)
            if swap_endian
                arr .= bswap.(arr)
            end
            return arr
        end

        # --- Read Header ---
        version = read_scalar(fid, Int32)
        if version < 20000
            @error "Error: this version of grid file not compatible. It only reads tortuga2 files."
            return nothing
        end

        nbrBlocks = read_scalar(fid, Int32)

        gridfile_bytes = read(fid, IO_NAME_LENGTH)
        gridfile = strip(String(gridfile_bytes), '\0')

        d = 3 # Dimension of space
        temp = Vector{Int32}(undef, d * nbrBlocks)
        read_array!(fid, temp)
        N = permutedims(reshape(temp, (d, nbrBlocks)))

        remainingHeaderSize = read_scalar(fid, Int32)
        if remainingHeaderSize > 0
            println("Skipping $remainingHeaderSize bytes of extra header")
            skip(fid, remainingHeaderSize)
        end

        # --- Read Data Fields ---
        varname_list = String[]
        vartype_list = Int[]
        fields = Any[]

        while !eof(fid)
            vtyp = read_scalar(fid, Int32)
            precision = read_scalar(fid, Int32)
            float_type = (precision == 4) ? Float32 : Float64

            vname_bytes = read(fid, IO_NAME_LENGTH)
            vname = strip(String(vname_bytes), '\0')
            
            push!(varname_list, vname)
            push!(vartype_list, vtyp)

            typeOfVariable = get(typekeyDict, vtyp, "UNKNOWN")
            println("--> Reading field $vname of type $typeOfVariable")

            local var 

            if vtyp == CV_SCALAR
                var = Vector{Array{float_type}}(undef, nbrBlocks)
                for blk in 1:nbrBlocks
                    Nx, Ny, Nz = N[blk, :]
                    temp = read_array!(fid, Vector{float_type}(undef, Nx * Ny * Nz))
                    var[blk] = permutedims(reshape(temp, (Nz, Ny, Nx)), (3, 2, 1))
                end

            elseif vtyp == CV_VECTOR
                var = Vector{Array{float_type}}(undef, nbrBlocks)
                for blk in 1:nbrBlocks
                    var[blk] = zeros(float_type, d, N[blk, :]...)
                end
                for l in 1:d, blk in 1:nbrBlocks
                    Nx, Ny, Nz = N[blk, :]
                    temp = read_array!(fid, Vector{float_type}(undef, Nx * Ny * Nz))
                    var[blk][l, :, :, :] = permutedims(reshape(temp, (Nz, Ny, Nx)), (3, 2, 1))
                end

            elseif vtyp == CV_TENSOR
                var = Vector{Array{float_type}}(undef, nbrBlocks)
                for blk in 1:nbrBlocks
                    var[blk] = zeros(float_type, d, d, N[blk, :]...)
                end
                for l1 in 1:d, l2 in 1:d, blk in 1:nbrBlocks
                    Nx, Ny, Nz = N[blk, :]
                    temp = read_array!(fid, Vector{float_type}(undef, Nx * Ny * Nz))
                    var[blk][l1, l2, :, :, :] = permutedims(reshape(temp, (Nz, Ny, Nx)), (3, 2, 1))
                end

            elseif vtyp == FA_SCALAR
                var = [Dict{String, Array{float_type}}() for _ in 1:nbrBlocks]
                faList = ["i", "j", "k"]
                for (faceDir, faKey) in enumerate(faList)
                    for blk in 1:nbrBlocks
                        Nb_vec = copy(N[blk, :])
                        Nb_vec[faceDir] += 1
                        Nx, Ny, Nz = Nb_vec
                        temp = read_array!(fid, Vector{float_type}(undef, Nx * Ny * Nz))
                        var[blk][faKey] = permutedims(reshape(temp, (Nz, Ny, Nx)), (3, 2, 1))
                    end
                end

            elseif vtyp == FA_VECTOR
                var = [Dict{String, Array{float_type}}() for _ in 1:nbrBlocks]
                faList = ["i", "j", "k"]
                for (faceDir, faKey) in enumerate(faList), blk in 1:nbrBlocks
                    Nb_vec = copy(N[blk, :])
                    Nb_vec[faceDir] += 1
                    var[blk][faKey] = zeros(float_type, d, Nb_vec...)
                end
                for l in 1:d, (faceDir, faKey) in enumerate(faList), blk in 1:nbrBlocks
                    Nb_vec = copy(N[blk, :])
                    Nb_vec[faceDir] += 1
                    Nx, Ny, Nz = Nb_vec
                    temp = read_array!(fid, Vector{float_type}(undef, Nx * Ny * Nz))
                    var[blk][faKey][l, :, :, :] = permutedims(reshape(temp, (Nz, Ny, Nx)), (3, 2, 1))
                end

            elseif vtyp == FA_TENSOR
                var = [Dict{String, Array{float_type}}() for _ in 1:nbrBlocks]
                faList = ["i", "j", "k"]
                for (faceDir, faKey) in enumerate(faList), blk in 1:nbrBlocks
                    Nb_vec = copy(N[blk, :])
                    Nb_vec[faceDir] += 1
                    var[blk][faKey] = zeros(float_type, d, d, Nb_vec...)
                end
                for l1 in 1:d, l2 in 1:d, (faceDir, faKey) in enumerate(faList), blk in 1:nbrBlocks
                    Nb_vec = copy(N[blk, :])
                    Nb_vec[faceDir] += 1
                    Nx, Ny, Nz = Nb_vec
                    temp = read_array!(fid, Vector{float_type}(undef, Nx * Ny * Nz))
                    var[blk][faKey][l1, l2, :, :, :] = permutedims(reshape(temp, (Nz, Ny, Nx)), (3, 2, 1))
                end

            else
                @error "Error: Cannot yet read variable of type $vtyp"
                break
            end
            push!(fields, var)
        end

        # --- Finalize and Return ---
        data = Dict{String, Any}(zip(varname_list, fields))
        datatype = Dict{String, Int}(zip(varname_list, vartype_list))

        return data, datatype, String(gridfile)
    end 
end