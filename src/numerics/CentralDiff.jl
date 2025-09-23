function CentralDiff(f, x)
    fx = similar(x)
    fx[2:end-1] = (f.(x[3:end]) - f.(x[1:end-2])) ./ (x[3:end] - x[1:end-2])
    fx[1] = (f(x[2]) - f(x[1])) / (x[2] - x[1])               
    fx[end] = (f(x[end]) - f(x[end-1])) / (x[end] - x[end-1]) 
    return fx
end