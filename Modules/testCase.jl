include("functionalPanel.jl")
# include("Modules/aeroModules.jl")

using .AeroModules
using BenchmarkTools

N_panels = 60
coords = open("S1223.dat") do f 
    [ Tuple(parse.(Float64,(split(line)))) for line in eachline(f) ] 
    end
x_ends, y_ends = AeroModules.cosineAirfoil([ x for (x, y) in coords ], [ y for (x, y) in coords ], N_panels)

# cosine = open("cosinepanels.dat", "w") do f
#     [ write(f, "$x $y \n") for (x, y) in zip(x_ends, y_ends) ]
# end

# Generate NACA airfoil coordinates
# N_panels = 100
# xs, ys = AeroModules.NACA4((2,4,1,2), 1.0, N_panels, true)
# coords = zip(xs, ys)
# x_ends, y_ends = AeroModules.cosineAirfoil([ x for (x, y) in coords ], [ y for (x, y) in coords ], N_panels)

# File coordinates
# x_ends, y_ends = [ x for (x, y) in coords ], [ y for (x, y) in coords ]

x_ends, y_ends = reverse(x_ends), reverse(y_ends)

panels = [ AeroModules.DoubletPanel2D((xs, ys), (xe, ye)) for (xs, xe, ys, ye) in zip(x_ends[1:end-1], x_ends[2:end], y_ends[1:end-1], y_ends[2:end]) ];

# Diamond wedge test
# panels = [ AeroModules.DoubletPanel2D((1.0, 0.0), (0.5, -0.1)), AeroModules.DoubletPanel2D((0.5, -0.1), (0.0, 0.0)), AeroModules.DoubletPanel2D((0.0, 0.0), (0.5, 0.1)), AeroModules.DoubletPanel2D((0.5, 0.1), (1.0, 0.0)) ]

# Hexagon
# panels = [ AeroModules.DoubletPanel2D((1.0, 0.0), (0.75, -0.1)), AeroModules.DoubletPanel2D((0.75, -0.1), (0.25, -0.1)), AeroModules.DoubletPanel2D((0.25, -0.1), (0.0, 0.0)), AeroModules.DoubletPanel2D((0.0, 0.0), (0.25, 0.1)), AeroModules.DoubletPanel2D((0.25, 0.1), (0.75, 0.1)), AeroModules.DoubletPanel2D((0.75, 0.1), (1.0, 0.0))]

uniform = AeroModules.Uniform2D(5.0, 5.0)

@time cl, kcl, cps = AeroModules.aeroCoefficients(panels, uniform)
println("Lift Coefficient: $cl, $kcl")
# println("Pressure Coefficients: $cps")