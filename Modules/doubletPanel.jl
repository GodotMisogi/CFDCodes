include("./doubletSourcePanelSolver.jl")
using .panelSolver
using BenchmarkTools
# using PyPlot
# using PyCall
# jtplot = pyimport("jupyterthemes.jtplot")
# jtplot.style(grid=false)
# rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
# rcParams["font.size"] = 16
# rcParams["font.family"] = "serif"
# rcParams["text.usetex"] = true

# Grid parameters
# x_domain, y_domain = (-1, 2), (-1, 1)
# grid_size = 50
# x_dom, y_dom = range(x_domain[1], length=grid_size, stop=x_domain[2]), range(y_domain[1], length=grid_size, stop=y_domain[2])
# X = repeat(x_dom', grid_size)
# Y = repeat(y_dom, 1, grid_size);

# Import airfoil coordinates file
N_panels = 60
filename = "S1223.dat"
coords = open(filename) do f 
    [ Tuple(parse.(Float64,(split(line)))) for line in eachline(f) ] 
    end
x_ends, y_ends = panelSolver.cosineAirfoil([ x for (x, y) in coords ], [ y for (x, y) in coords ], N_panels)

# Generate NACA airfoil coordinates
# N_panels = 160
# xs, ys = panelSolver.NACA4((4,4,1,2), 1.0, N_panels, true)
# coords = zip(xs, ys)
# x_ends, y_ends = panelSolver.cosineAirfoil([ x for (x, y) in coords ], [ y for (x, y) in coords ], N_panels)

# Generate flat plate coordinates
# coords = zip(reverse(range(0, 1, length=10)), zeros(10));
# x_ends, y_ends = panelSolver.cosinePanels([ x for (x, y) in coords ], [ y for (x, y) in coords], 30);

# Save cosine-distribution interpolated points
cosine = open("cosinepanels.dat", "w") do f
    [ write(f, "$x $y \n") for (x, y) in zip(x_ends, y_ends) ]
end

# File coordinates
# x_ends, y_ends = [ x for (x, y) in coords ], [ y for (x, y) in coords ]

# Reverse coordinates to follow clockwise convention
x_ends, y_ends = reverse(x_ends), reverse(y_ends)

# Using airfoil panels
panels = [ panelSolver.DoubletSourcePanel2D(xs, ys, xe, ye) for (xs, xe, ys, ye) in zip(x_ends[1:end-1], x_ends[2:end], y_ends[1:end-1], y_ends[2:end]) ];

# Diamond wedge test case
# panels = [ panelSolver.DoubletSourcePanel2D(1.0, 0.0, 0.5, -0.1), panelSolver.DoubletSourcePanel2D(0.5, -0.1, 0.0, 0.0), panelSolver.DoubletSourcePanel2D(0.0, 0.0, 0.5, 0.1), panelSolver.DoubletSourcePanel2D(0.5, 0.1, 1.0, 0.0) ]

# # Hexagon
# panels = [ panelSolver.DoubletSourcePanel2D(1.0, 0.0, 0.75, -0.1), panelSolver.DoubletSourcePanel2D(0.75, -0.1, 0.25, -0.1), panelSolver.DoubletSourcePanel2D(0.25, -0.1, 0.0, 0.0), panelSolver.DoubletSourcePanel2D(0.0, 0.0, 0.25, 0.1), panelSolver.DoubletSourcePanel2D(0.25, 0.1, 0.75, 0.1), panelSolver.DoubletSourcePanel2D(0.75, 0.1, 1.0, 0.0) ]

# Uniform2D construction
uniform_mag, uniform_ang = 5.0, 5.0
uniform = panelSolver.Uniform2D(uniform_mag, uniform_ang);


# Print panel attributes
# [ print("Start: ", (panel.xs, panel.ys), ", End: ", (panel.xe, panel.ye), '\n', "Normal: ", panel.normal, ", Tangent: ", panel.tangent, ", Angle: ", rad2deg(panel.angle), ", Position: ", panel.loc, '\n') for panel in panels]
# print(panels[1].length, " ", panels[end].length, '\n')

# Solve system
@time airfoil, strengths, cps, cl, error = panelSolver.panelSolver2D(panels, uniform, false, true);

# For plotting purposes
# vels, pots = panelSolver.gridData(airfoil, zip(X,Y));

# Pressure coefficients
# cp = panelSolver.pressureCoefficient2D(vels, uniform.magnitude);
# print([ panel.source_strength for panel in panels ])
println("Lift Coefficient: ", cl)
# DEBUGGING
# println("Strengths:", strengths)
# println("Pressure Coefficients:", cps)
# println("'Error': ", error)

# Airfoil plotter
# fig1 = figure(1, dpi=300)
# xlabel("\$x\$")
# ylabel("\$y\$")
# plot([ x for (x,y) in coords], [ y for (x,y) in coords ], linestyle="--", linewidth=1)
# plot([ panel.xs for panel in panels ], [ panel.ys for panel in panels ], linestyle="-", linewidth=0.5, marker=".", markersize=3, color="#CD2305")
# axis("scaled", adjustable="box")
# # xlim(-0.02,0.02)
# ylim(minimum([ y for (x,y) in coords ]) - 0.05, maximum([ y for (x,y) in coords ]) + 0.05)
# grid()

# Plot pressure coefficient
# fig2 = figure(2, figsize=(10, 5), dpi=300)
# grid()
# plot([p.xc for p in panels if p.loc == "upper" ], [p.cp for p in panels if p.loc == "upper" ], marker=".", label="Upper")
# plot([p.xc for p in panels if p.loc == "lower" ], [p.cp for p in panels if p.loc == "lower" ], marker=".", label="Lower")
# xlabel("\$x\$")
# ylabel("\$C_{p_{local}}\$")
# xlim(0, 1.0)
# ylim(-2,2)
# ylim(maximum([p.cp for p in panels]) + 0.05, minimum([p.cp for p in panels]) - 0.05)
# tight_layout()
# legend()
# show()

# Plot flowfield
# width = 8.0
# height = (y_domain[2] - y_domain[1]) / (x_domain[2] - x_domain[1]) * width
# fig3 = figure(3, figsize=(1.1width, height), dpi=300)
# streamplot(X, Y, first.(vels), last.(vels), density=2)
# # scatter([p.xc for p in source_panels], [p.yc for p in source_panels], marker=".", label="Center-points", color="white", zorder=3)
# contourf(X, Y, cp, 50, cmap=PyPlot.cm.rainbow)
# # colorbar(label="\$\\phi_{\\mathrm{panel}}\$")
# colorbar(label="\$C_p\$")
# # fill([x for (x,y) in coords], [y for (x,y) in coords], color="k", zorder=3)
# plot([p.xs for p in panels], [p.ys for p in panels], linestyle="--", linewidth=0.5, label="Endpoints", color="orange", zorder=3)
# xlim(x_domain)
# ylim(y_domain)
# title("Number of Doublet-Source Panels: {$N_panels}")
# tight_layout();
