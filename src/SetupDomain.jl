function SetupDomain(airfoil, radius, vertN, horzN; type = "cgrid")
    # add a point to be cut off 
    horzN = horzN + 1  
    # throw away the last column
    airfoil = airfoil[:, 1:2]

    if type == "cgrid"
        # build the C-grid around the airfoil
        airfoilBoundary = BuildCGrid(airfoil, radius, vertN)
        
        # build out the branch cut on the left side
        branchCutLeft = BuildBranchCut(airfoilBoundary, radius, vertN, horzN; direction="left" )
        
        # build out the branch cut on the right side
        branchCutRight = BuildBranchCut(airfoilBoundary, radius, vertN, horzN; direction="right")

        # combine the airfoil boundary and the branch cuts, removing the duplicate points
        boundaryTop = vcat(branchCutLeft[1][1:end-1,:], airfoilBoundary[1], branchCutRight[1][2:end,:])

        boundaryBottom = vcat(branchCutLeft[3][1:end-1,:], airfoilBoundary[3], branchCutRight[3][2:end,:])

        boundaryLeft = branchCutLeft[4]
        boundaryRight = branchCutRight[2]
    else 
        error("Invalid type. Use 'cgrid'.")
    end

    return [boundaryTop, boundaryRight, boundaryBottom, boundaryLeft]
end



function BuildBranchCut(airfoilBoundary, radius, vertN, horzN; direction)
    #---------------------------
    # Description: Build the branch cut from the trailing edge of the airfoil to the outer boundary of size radius.
    # Input: airfoilBoundary, radius
    # Output: branchCut [top, right, bottom, left] 
    #---------------------------
    # Set the number of points for the branch cut

    # get the top and bottom edges of the airfoil boundary
    top = airfoilBoundary[1]
    bottom = airfoilBoundary[3]


        
    tV = range(0, 1; length=vertN)
    tH = range(0, 1; length=horzN)

    trailingEdge = bottom[1,:]
    xStartBottom, yStartBottom = [radius, trailingEdge[2]]
    xEndBottom, yEndBottom = trailingEdge

    branchCutBottom = hcat(
        xStartBottom .+ tH .* (xEndBottom - xStartBottom),   # x‐coords
        yStartBottom .+ tH .* (yEndBottom - yStartBottom)    # y‐coords
    )
    
    if direction == "left"
        # build a straight line from (radius, 0) to (trailingEdge[1], trailingEdge[2])
        # build the straight line, nPoints×2 array

        # get the position of the top edge
        topEdge = top[1,:]
        # build a straight line from (radius, -radius) to (topEdge[1], topEdge[2])
        xStartTop, yStartTop = [radius, -radius]
        xEndTop, yEndTop = topEdge
    elseif direction == "right"
        # trailingEdge = bottom[end,:]
        # # build a straight line from (trailingEdge[1], trailingEdge[2]) to (radius, 0)
        # xStartBottom, yStartBottom = trailingEdge
        # xEndBottom, yEndBottom = [radius, trailingEdge[2]]

        # reverse the bottom edge
        branchCutBottom = reverse(branchCutBottom, dims=1)
        
        # get the position of the top edge
        topEdge = top[end,:]
        # build a straight line from (topEdge[1], topEdge[2]) to (radius, radius)
        xStartTop, yStartTop = topEdge
        xEndTop, yEndTop = [radius, radius]

    else
        error("Invalid direction. Use 'left' or 'right'.")
    end


    

    # build the straight line, nPoints×2 array
    branchCutTop = hcat(
      xStartTop .+ tH .* (xEndTop - xStartTop),   # x‐coords
      yStartTop .+ tH .* (yEndTop - yStartTop)    # y‐coords
    )
    
    if direction == "left"
        # connect the top and bottom edges on the left side
        xStartLeft, yStartLeft = branchCutBottom[1, :]
        xEndLeft, yEndLeft = branchCutTop[1, :]
        branchCutLeft = hcat(
            xStartLeft .+ tV .* (xEndLeft - xStartLeft),   # x‐coords
            yStartLeft .+ tV .* (yEndLeft - yStartLeft)    # y‐coords
        )


        # use the airfoil boundary to get the right edge
        branchCutRight = airfoilBoundary[4]
    elseif direction == "right"
        # use the airfoil boundary to get the left edge
        branchCutLeft = airfoilBoundary[2]

        # connect the top and bottom edges on the right side
        xStartRight, yStartRight = branchCutBottom[end, :]
        xEndRight, yEndRight = branchCutTop[end, :]
        branchCutRight = hcat(
            xStartRight .+ tV .* (xEndRight - xStartRight),   # x‐coords
            yStartRight .+ tV .* (yEndRight - yStartRight)    # y‐coords
        )
    else 
        error("Invalid direction. Use 'left' or 'right'.")
    end

    # return the branch cut as a 4-element array
    return [branchCutTop, branchCutRight, branchCutBottom, branchCutLeft]
end

function BuildCGrid(airfoil, radius, vertN)
    #---------------------------
    # Description: Extend the boundary around the airfoil to radius size using the normal vector. The function calls ClampToDomain! to ensure that the outer domain is clamped to the desired radius. 
    # Input: airfoil, radius 
    # Output: clamped outer domain [top, right, bottom, left]
    #---------------------------
    airfoilSize = size(airfoil, 1)
    # allocate
    tangent = zeros(airfoilSize, 2)
    normal  = zeros(airfoilSize,  2)

    # compute central differences along the airfoil
    dx = airfoil[3:end,1] .- airfoil[1:end-2,1]
    dy = airfoil[3:end,2] .- airfoil[1:end-2,2]
    
    # save in tangent vector
    tangent[2:end-1,1] = dx
    tangent[2:end-1,2] = dy

    # one sided differences for the first and last point
    tangent[1, :] = airfoil[2, :] .- airfoil[1, :]
    tangent[end, :] = airfoil[end, :] .- airfoil[end-1, :]

    # normalize the tangent vector
    tangent = tangent ./ sqrt.(tangent[:,1].^2 .+ tangent[:,2].^2)

    # compute the normal vector
    normal[:,1] = -tangent[:,2]
    normal[:,2] = tangent[:,1]

    offsets = radius .* normal

    # add back to the airfoil points
    outerDomain = airfoil .+ offsets
    
    # clamp the outer domain to the desired cgrid radius. Adding a straight away
    ClampToDomain!(outerDomain, radius)
    
    # save the bottom and top edges
    top = outerDomain
    bottom = airfoil

    # connect top and bottom with a line
    left = zeros(vertN, 2)
    x_start, y_start = bottom[1, :]
    x_end, y_end = top[1, :]
    
    t = range(0, 1; length=vertN)
    
    # build the straight line, vert_N×2 array
    left = hcat(
      x_start .+ t .* (x_end - x_start),   # x‐coords
      y_start .+ t .* (y_end - y_start)    # y‐coords
    )
    
    right = zeros(vertN, 2)
    x_start, y_start = bottom[end, :]
    x_end, y_end = top[end, :]
    right = hcat(
      x_start .+ t .* (x_end - x_start),   # x‐coords
      y_start .+ t .* (y_end - y_start)    # y‐coords
    )
    return [top, right, bottom, left]
end

function ClampToDomain!(outerDomain, radius)
    # clamp the outer domain to the desired cgrid
    # loop through the points
    for i in eachindex(outerDomain[:, 1])
        x, y = outerDomain[i, :]
        # clamp to straight away under the airfoil
        if x > 0 && y < 0        
            y = -radius
            outerDomain[i, :] = [x, y]
        # Clamp to the straight away above the airfoil
        elseif x > 0 && y > 0
            y = radius
            outerDomain[i, :] = [x, y]
        # otherwise clamp to the semi circle
        else
            # compute original distance 
            d = sqrt(x^2 + y^2)
            # scale the vector to the radius
            x = x * (radius / d)
            y = y * (radius / d)

            outerDomain[i, :] = [x,y]
        end
    end
end