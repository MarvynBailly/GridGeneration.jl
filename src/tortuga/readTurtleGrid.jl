"""
    readTurtleGrid(gridfilename::String)

Reads grid data from a "Turtle" format file, correctly handling row-major
data layout for node coordinates. It uses dictionaries to store boundary and
interface information.

# Arguments
- `gridfilename::String`: Path to the grid file.

# Returns
- `Xno`, `Xcv`, `Xfa`, `interfaceInfo`, `bndInfo`
"""
function readTurtleGrid(gridfilename::String)
    # Parameters from Turtle_parameters.h
    MAGICNUMBER = 751123
    IO_NAME_LENGTH = 200

    println("Reading $gridfilename")
    if occursin("reduced", gridfilename) || occursin("wallgrid", gridfilename)
        println("NOTE: reduced or wallgrid files do not contain any boundary info and are for post-processing only")
    end

    open(gridfilename, "r") do fid
        # --- Endianness Check and Helper Functions ---
        magicno = read(fid, Int32)
        swap_endian = false
        if magicno != MAGICNUMBER
            magicno = bswap(magicno)
            if magicno == MAGICNUMBER
                swap_endian = true
                println("File is non-native endian. Swapping byte order.")
            else
                @error "Cannot find correct endian. Terminating."
                return
            end
        end

        read_scalar(io, type) = swap_endian ? bswap(read(io, type)) : read(io, type)

        function read_array!(io, arr::AbstractArray)
            read!(io, arr)
            if swap_endian
                arr .= bswap.(arr)
            end
        end

        # --- Read Header ---
        version = read_scalar(fid, Int32)
        if version < 20000
            @error "This version of grid file is not compatible. It only reads tortuga2 files."
            return
        end

        nbrBlocks, nbrInterfaces, nbrBoundaries = [read_scalar(fid, Int32) for _ in 1:3]
        println("--> The grid has $nbrBlocks blocks, $nbrInterfaces interfaces, and $nbrBoundaries boundaries")

        temp = Vector{Int32}(undef, 3 * nbrBlocks)
        read_array!(fid, temp)
        N = permutedims(reshape(temp, (3, nbrBlocks)))

        interfaceInfo = Any[]
        bndInfo = Any[]

        # --- Read Interface and Boundary Data based on Version ---
        if version < 30101 # Version < 3.1
            println("--> Reading interface and boundaries in block number-IJK format (< 3.1)")
            for _ in 1:nbrInterfaces
                interface_dict = Dict{String, Any}()
                interface_dict["blockA"] = read_scalar(fid, Int32)
                IJK_A = [read_scalar(fid, Int32) for _ in 1:3]
                b = Int(interface_dict["blockA"]) + 1
                sign_A = IJK_A .% 2
                dir_A = IJK_A .÷ 2 .+ 1
                interface_dict["start_blkA"] = zeros(Int32, 3)
                interface_dict["end_blkA"] = zeros(Int32, 3)
                sign_A_flipped = 1 .- sign_A
                interface_dict["end_blkA"][dir_A] = N[b, dir_A] .* sign_A
                interface_dict["start_blkA"][dir_A] = N[b, dir_A] .* sign_A_flipped
                interface_dict["blockB"] = read_scalar(fid, Int32)
                IJK_B = [read_scalar(fid, Int32) for _ in 1:3]
                if IJK_B[1] % 2 == 0; IJK_B[1] += 1; else; IJK_B[1] -= 1; end
                b = Int(interface_dict["blockB"]) + 1
                sign_B = IJK_B .% 2
                dir_B = IJK_B .÷ 2 .+ 1
                interface_dict["start_blkB"] = zeros(Int32, 3)
                interface_dict["end_blkB"] = zeros(Int32, 3)
                sign_B_flipped = 1 .- sign_B
                interface_dict["end_blkB"][dir_B] = N[b, dir_B] .* sign_B
                interface_dict["start_blkB"][dir_B] = N[b, dir_B] .* sign_B_flipped
                interface_dict["offset"] = [read_scalar(fid, Float64) for _ in 1:3]
                push!(interfaceInfo, interface_dict)
            end
            for _ in 1:nbrBoundaries
                bnd_dict = Dict{String, Any}()
                name_bytes = read(fid, IO_NAME_LENGTH)
                bnd_dict["name"] = strip(split(String(name_bytes), '\0')[1])
                bnd_dict["nbrFaces"] = read_scalar(fid, Int32)
                face_vector = Dict{String, Any}[]
                for _ in 1:bnd_dict["nbrFaces"]
                    face_dict = Dict{String, Any}()
                    bse_data = [read_scalar(fid, Int32) for _ in 1:4]
                    face_dict["block"] = bse_data[1]
                    IJK = bse_data[2:4]
                    b = Int(face_dict["block"]) + 1
                    sign_val = IJK .% 2
                    dir = IJK .÷ 2 .+ 1
                    face_dict["start"] = zeros(Int32, 3)
                    face_dict["end"] = zeros(Int32, 3)
                    sign_flipped = 1 .- sign_val
                    face_dict["end"][dir] = N[b, dir] .* sign_val
                    face_dict["start"][dir] = N[b, dir] .* sign_flipped
                    push!(face_vector, face_dict)
                end
                bnd_dict["faceInfo"] = face_vector
                push!(bndInfo, bnd_dict)
            end
        else # Version >= 3.1
            is_newest_format = version >= 30147
            if is_newest_format; println("--> Reading in patch format (>= 3.1.47)"); else; println("--> Reading in patch format (>= 3.1)"); end
            for _ in 1:nbrInterfaces
                interface_dict = Dict{String, Any}()
                interface_dict["blockA"] = read_scalar(fid, Int32)
                start_end_dir_A = [read_scalar(fid, Int32) for _ in 1:6]
                interface_dict["blockB"] = read_scalar(fid, Int32)
                start_end_dir_B = [read_scalar(fid, Int32) for _ in 1:6]
                interface_dict["offset"] = [read_scalar(fid, Float64) for _ in 1:3]
                interface_dict["angle"] = is_newest_format ? read_scalar(fid, Float64) : 0.0
                interface_dict["start_blkA"] = start_end_dir_A[1:3]
                interface_dict["end_blkA"] = start_end_dir_A[4:6]
                interface_dict["start_blkB"] = start_end_dir_B[1:3]
                interface_dict["end_blkB"] = start_end_dir_B[4:6]
                push!(interfaceInfo, interface_dict)
            end
            for _ in 1:nbrBoundaries
                bnd_dict = Dict{String, Any}()
                name_bytes = read(fid, IO_NAME_LENGTH)
                bnd_dict["name"] = strip(split(String(name_bytes), '\0')[1])
                bnd_dict["nbrFaces"] = read_scalar(fid, Int32)
                face_vector = Dict{String, Any}[]
                for _ in 1:bnd_dict["nbrFaces"]
                    face_dict = Dict{String, Any}()
                    data = [read_scalar(fid, Int32) for _ in 1:7]
                    face_dict["block"] = data[1]
                    face_dict["start"] = data[2:4]
                    face_dict["end"] = data[5:7]
                    push!(face_vector, face_dict)
                end
                bnd_dict["faceInfo"] = face_vector
                push!(bndInfo, bnd_dict)
            end
        end

        # --- Read Node Coordinates ---
        d = 3
        Xno = [zeros(Float64, d, N[blk, 1] + 1, N[blk, 2] + 1, N[blk, 3] + 1) for blk in 1:nbrBlocks]
        
        # This now handles the row-major to column-major conversion correctly.
        for l in 1:d, blk in 1:nbrBlocks
            # Get the dimensions for this block
            Nx, Ny, Nz = N[blk, :] .+ 1

            # Read the flat data
            count = Nx * Ny * Nz
            temp = Vector{Float64}(undef, count)
            read_array!(fid, temp)

            # 1. Reshape with reversed dimensions to match memory layout
            temp_reshaped = reshape(temp, (Nz, Ny, Nx))

            # 2. Permute dimensions to get the correct logical order (d, Nx, Ny, Nz)
            Xno[blk][l, :, :, :] = permutedims(temp_reshaped, (3, 2, 1))
        end

        if !eof(fid); @warn "End of file not reached after reading all expected data."; end

        # --- Compute Face and Cell Centers ---
        Xcv = [zeros(Float64, d, N[blk, 1], N[blk, 2], N[blk, 3]) for blk in 1:nbrBlocks]
        Xfa = [Dict{String, Array{Float64}}() for _ in 1:nbrBlocks]

        for blk in 1:nbrBlocks
            @views begin
                Xno_blk = Xno[blk]
                Xcv[blk] = 0.125 * (Xno_blk[:, 1:end-1, 1:end-1, 1:end-1] .+ Xno_blk[:, 2:end, 1:end-1, 1:end-1] .+ Xno_blk[:, 1:end-1, 2:end, 1:end-1] .+ Xno_blk[:, 1:end-1, 1:end-1, 2:end] .+ Xno_blk[:, 2:end, 2:end, 1:end-1] .+ Xno_blk[:, 2:end, 1:end-1, 2:end] .+ Xno_blk[:, 1:end-1, 2:end, 2:end] .+ Xno_blk[:, 2:end, 2:end, 2:end])
                Xfa[blk]["i"] = 0.25 * (Xno_blk[:, :, 1:end-1, 1:end-1] .+ Xno_blk[:, :, 1:end-1, 2:end] .+ Xno_blk[:, :, 2:end, 1:end-1] .+ Xno_blk[:, :, 2:end, 2:end])
                Xfa[blk]["j"] = 0.25 * (Xno_blk[:, 1:end-1, :, 1:end-1] .+ Xno_blk[:, 1:end-1, :, 2:end] .+ Xno_blk[:, 2:end, :, 1:end-1] .+ Xno_blk[:, 2:end, :, 2:end])
                Xfa[blk]["k"] = 0.25 * (Xno_blk[:, 1:end-1, 1:end-1, :] .+ Xno_blk[:, 1:end-1, 2:end, :] .+ Xno_blk[:, 2:end, 1:end-1, :] .+ Xno_blk[:, 2:end, 2:end, :])
            end
        end

        return Xno, Xcv, Xfa, interfaceInfo, bndInfo
    end
end



"""
    ImportTurtleGrid(gridfilename::String)

A simplified interface to read a Turtle grid file and extract 2D node coordinates
along with interface and boundary information.

# Arguments
- `gridfilename::String`: Path to the grid file.

# Returns
- `Xno2D`, `interfaceInfo`, `bndInfo`
""" 
function ImportTurtleGrid(gridfilename::String)
    Xno, Xcv, Xfa,interfaceInfo,bndInfo = readTurtleGrid(gridfilename)
    
    Xno2D = [Xno_blk[1:2,:,:,1] for Xno_blk in Xno]
    Xcv2D = [Xcv_blk[1:2,:,:,1] for Xcv_blk in Xcv]

    return Xno2D, Xcv2D, Xfa, interfaceInfo, bndInfo
end