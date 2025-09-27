function GetRectangleDomain()
    # Define height and width of the rectangle
    height = 2.0
    width  = 4.0

    # Define number of points
    num_points_height = 500
    num_points_width  = 500

    top = hcat(range(0, stop=width, length=num_points_width), fill(height, num_points_width))
    bottom = hcat(range(0, stop=width, length=num_points_width), fill(0.0, num_points_width))
    right = hcat(fill(width, num_points_height), range(0, stop=height, length=num_points_height))
    left = hcat(fill(0.0, num_points_height), range(0, stop=height, length=num_points_height))

    return [top, right, bottom, left]
end