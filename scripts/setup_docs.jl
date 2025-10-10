import Pkg
pkgpath = abspath(joinpath(@__DIR__, ".."))
println("Developing local package from: ", pkgpath)
Pkg.develop(Pkg.PackageSpec(path=pkgpath))
Pkg.instantiate()
Pkg.precompile()
println("Docs project setup complete")
