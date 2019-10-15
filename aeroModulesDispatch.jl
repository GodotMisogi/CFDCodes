# Implementing panel methods using multiple dispatch. Warning: Slow AF

module panelSolver

    using QuadGK

    abstract type Solution end

    struct Source2D <: Solution
        strength :: Float64
        x0 :: Float64
        y0 :: Float64
        function Source2D(strength, x0, y0)

        end
    end

    velocity(source :: Source2D, x :: Float64, y :: Float64) = (source.strength/(2π)*(x - source.x0)/((x - source.x0)^2 + (y - source.y0)^2), source.strength/(2π)*(y - source.y0)/((x - source.x0)^2 + (y - source.y0)^2))
    stream(source :: Source2D, x :: Float64, y :: Float64) = source.strength/(2π)*atan(y - source.y0, x - source.x0)
    potential(source :: Source2D, x :: Float64, y :: Float64) = source.strength/(4π)*log((x - source.x0)^2 + (y - source.y0)^2)

    struct Uniform2D <: Solution
        magnitude :: Float64
        angle :: Float64
        deg :: Float64
        function Uniform2D(magnitude, deg)
            angle = deg2rad(deg)
            new(magnitude, angle, deg)
        end
    end
    
    velocity(uniform :: Uniform2D, x :: Float64, y :: Float64) = (uniform.magnitude*cos(uniform.angle), uniform.magnitude*sin(uniform.angle))
    stream(uniform :: Uniform2D, x :: Float64, y :: Float64) = uniform.magnitude*(y*cos(uniform.angle) - x*sin(uniform.angle))
    potential(uniform :: Uniform2D, x :: Float64, y :: Float64) = uniform.magnitude*(x*cos(uniform.angle) + y*sin(uniform.angle))

    struct Doublet2D <: Solution
        strength :: Float64
        x0 :: Float64
        y0 :: Float64
    end

    velocity(doublet :: Doublet2D, x :: Float64, y :: Float64) = (-doublet.strength/(2π)*((x - doublet.x0)^2 - (y - doublet.y0)^2)/((x - doublet.x0)^2 + (y - doublet.y0)^2)^2, -doublet.strength/(2π)*2*(x - doublet.x0)*(y - doublet.y0)/((x - doublet.x0)^2 + (y - doublet.y0)^2)^2)
    stream(doublet :: Doublet2D, x :: Float64, y :: Float64) = -doublet.strength/(2π)*(y - doublet.y0)/((x - doublet.x0)^2 + (y - doublet.y0)^2)
    potential(doublet :: Doublet2D, x :: Float64, y :: Float64) = -doublet.strength/(2π)*(x - doublet.x0)/((x - doublet.x0)^2 + (y - doublet.y0)^2)

    struct Vortex2D <: Solution
        strength :: Float64
        x0 :: Float64
        y0 :: Float64
    end

    velocity(vortex :: Vortex2D, x :: Float64, y :: Float64) = (-vortex.strength/(2π)*(y - vortex.y0)/((x - vortex.x0)^2 + (y - vortex.y0)^2), vortex.strength/(2π)*(x - vortex.x0)/((x - vortex.x0)^2 + (y - vortex.y0)^2))
    stream(vortex :: Vortex2D, x :: Float64, y :: Float64) = -vortex.strength/(4π)*log((x - vortex.x0)^2 + (y - vortex.y0)^2)
    potential(vortex :: Vortex2D, x :: Float64, y :: Float64) = vortex.strength/(2π)*atan(y - vortex.y0, x - vortex.x0)

    abstract type Panel <: Solution end

    mutable struct SourcePanel2D <: Panel
        xs :: Float64
        ys :: Float64
        xe :: Float64 
        ye :: Float64
        xc :: Float64
        yc :: Float64
        length :: Float64
        angle :: Float64
        loc :: T where T <: AbstractString
        strength :: Float64
        vt :: Float64
        cp :: Float64
        function SourcePanel2D(xs, ys, xe, ye, strength=0.0, vt=0.0, cp=0.0)
            xc, yc   = (xe + xs)/2.0, (ye + ys)/2.0
            length = mag([ xe - xs, ye - ys ])
            angle = xe <= xs ? acos((ye - ys)/length) : π + acos(-(ye - ys)/length)
            loc = angle <= π ? "upper" : "lower"
            new(xs, ys, xe, ye, xc, yc, length, angle, loc, strength, vt, cp)
        end 
    end

    mutable struct VortexSourcePanel2D <: Panel
        xs :: Float64
        ys :: Float64
        xe :: Float64 
        ye :: Float64
        xc :: Float64
        yc :: Float64
        length :: Float64
        angle :: Float64
        loc :: T where T <: AbstractString
        strength :: Float64
        vt :: Float64
        cp :: Float64
        function VortexSourcePanel2D(xs, ys, xe, ye, strength=0.0, vt=0.0, cp=0.0)
            xc, yc   = (xe + xs)/2.0, (ye + ys)/2.0
            length = mag([ xe - xs, ye - ys ])
            angle = xe <= xs ? acos((ye - ys)/length) : π + acos(-(ye - ys)/length)
            loc = angle <= π ? "upper" : "lower"
            new(xs, ys, xe, ye, xc, yc, length, angle, loc, strength, vt, cp)
        end 
    end

    struct Panel2D <: Panel
        xs :: Float64
        ys :: Float64
        xe :: Float64 
        ye :: Float64
        xc :: Float64
        yc :: Float64
        length :: Float64
        angle :: Float64
        loc :: T where T <: AbstractString
        function Panel2D(xs, ys, xe, ye)
            xc, yc   = (xe + xs)/2.0, (ye + ys)/2.0
            length = mag([ xe - xs, ye - ys ])
            angle = xe <= xs ? acos((ye - ys)/length) : π + acos(-(ye - ys)/length)
            loc = angle <= π ? "upper" : "lower"
            new(xs, ys, xe, ye, xc, yc, length, angle, loc)
        end 
    end

    integral(panel :: Panel, x :: Float64, y :: Float64, dxdz :: Float64, dydz :: Float64) = quadgk(s -> ((x - (panel.xs - s*sin(panel.angle)))*dxdz + (y - (panel.ys + s*cos(panel.angle)))*dydz)/((x - (panel.xs - s*sin(panel.angle)))^2 + (y - (panel.ys + s*cos(panel.angle)))^2), 0., panel.length, atol=1e-10)[1]

    potential(panel :: Panel, x :: Float64, y :: Float64) = quadgk( (s) -> log((x - (panel.xs - s*sin(panel.angle)))^2 + (y - (panel.ys + s*cos(panel.angle)))^2), 0, panel.length, atol=1e-10)[1]
    
    velocity(panel :: Panel, x :: Float64, y :: Float64) = panel.strength/(2π).*(integral(panel, x, y, 1., 0.), integral(panel, x, y, 0., 1.))

    function sourcePanelSolver2D(panels :: Array{SourcePanel2D, 1}, uniform :: Uniform2D)
        """
        Applies the source panel method to a list of panels and a uniform flow.
        """
        strengths = solveStrengths(panels, uniform)
        vts = tangentialVelocities(panels, uniform)
        error = errorCalc(panels)

        return strengths, vts, error
    end

    function solveStrengths(panels :: Array{SourcePanel2D, 1}, uniform :: Uniform2D)
        """
        Solves for the source strengths of all panels.
        """
        # Construct matrix for normal direction.
        An = [ i == j ? 0.5 : 0.5/π*integral(panel_i, panel_j.xc, panel_j.yc, cos(panel_j.angle), sin(panel_j.angle)) for (j, panel_j) in enumerate(panels), (i, panel_i) in enumerate(panels) ]

        # Construct freestream RHS.
        bn = @. -uniform.magnitude*cos([ uniform.angle - panel.angle for panel in panels ]) 
        
        # Solve system.
        sources = An\bn

        # Update panel strengths.
        for (panel, strength) in zip(panels, sources)
            panel.strength = strength
        end
        return sources
    end
    
    function tangentialVelocities(panels :: Array{SourcePanel2D, 1}, uniform :: Uniform2D, pressure = true)
        """
        Solves for the velocities of all panels and their pressure coefficients.
        """
        # Construct matrix for tangential direction.
        At = [ i == j ? 0. : 0.5/π*integral(panel_i, panel_j.xc, panel_j.yc, -sin(panel_j.angle), cos(panel_j.angle)) for (j, panel_j) in enumerate(panels), (i, panel_i) in enumerate(panels) ]

        # Construct freestream RHS.
        bt = uniform.magnitude*sin.([ uniform.angle - panel.angle for panel in panels ]) 

        # Solve system.
        vts = At*[ panel.strength for panel in panels ] .+ bt

        # Update panel velocities and pressure coefficients.
        for (panel, vt) in zip(panels, vts)
            panel.vt = vt
            panel.cp = pressure ? pressureCoefficient2D(0., vt, uniform.magnitude) : nothing
        end

        return vts
    end
   
    function vortexSourcePanelSolver2D(panels :: Array{VortexSourcePanel2D, 1}, uniform :: Uniform2D)
        """
        Applies the vortex source panel method to a list of panels and a uniform flow.
        """
        strengths = solveStrengths(panels, uniform)
        vts = tangentialVelocities(strengths, panels, uniform)
        cl = liftCoefficient(strengths[end], panels, uniform)
        error = errorCalc(panels)

        return strengths, vts, cl, error
    end

    function solveStrengths(panels :: Array{VortexSourcePanel2D, 1}, uniform :: Uniform2D)
        """
        Solves for the source strengths of all panels.
        """
        num_panels = length(panels)

        # Construct source matrix for normal direction.
        source_matrix = [ i == j ? 0.5 : 0.5/π*integral(panel_i, panel_j.xc, panel_j.yc, cos(panel_j.angle), sin(panel_j.angle)) for (j, panel_j) in enumerate(panels), (i, panel_i) in enumerate(panels) ]

        # Vortex matrix
        vortex_matrix = [ i == j ? 0.0 : -0.5/π*integral(panel_i, panel_j.xc, panel_j.yc, sin(panel_j.angle), -cos(panel_j.angle)) for (j, panel_j) in enumerate(panels), (i, panel_i) in enumerate(panels) ]

        # Kutta condition
        kutta = zeros(num_panels + 1)
        kutta[1:end-1] = vortex_matrix[1, :] + vortex_matrix[end, :]
        kutta[end] = -sum(source_matrix[1, :] + source_matrix[end, :]) 
        
        # Global matrix construction
        An = zeros(num_panels + 1, num_panels + 1)
        An[1:end-1, 1:end-1] = source_matrix
        An[1:end-1, end] = sum(vortex_matrix, dims=2)
        An[end, :] = kutta

        # Construct freestream RHS.
        bn = push!(-uniform.magnitude .* cos.([ uniform.angle - panel.angle for panel in panels ]), -uniform.magnitude * (sin(uniform.angle - panels[1].angle) + sin(uniform.angle - panels[end].angle)))

        # Solve system.
        strengths = An\bn

        # Update panel strengths.
        for (panel, strength) in zip(panels, strengths)
            panel.strength = strength
        end

        return strengths
    end

    function tangentialVelocities(strengths :: Array{<: Real, 1}, panels :: Array{VortexSourcePanel2D, 1}, uniform :: Uniform2D, pressure = true)
        """
        Solves for the velocities of all panels and their pressure coefficients.
        """
        num_panels = length(panels)

        # Construct source matrix for normal direction.
        source_matrix = [ i == j ? 0.5 : 0.5/π*integral(panel_i, panel_j.xc, panel_j.yc, cos(panel_j.angle), sin(panel_j.angle)) for (j, panel_j) in enumerate(panels), (i, panel_i) in enumerate(panels) ]

        # Vortex matrix
        vortex_matrix = [ i == j ? 0.0 : -0.5/π*integral(panel_i, panel_j.xc, panel_j.yc, sin(panel_j.angle), -cos(panel_j.angle)) for (j, panel_j) in enumerate(panels), (i, panel_i) in enumerate(panels) ]

        At = zeros(num_panels, num_panels + 1)

        # Construct matrix for tangential direction.
        At[:, 1:end-1] = vortex_matrix
        At[:, end] = -sum(source_matrix, dims=2)
        
        # Construct freestream RHS.
        bt = @. uniform.magnitude*sin([ uniform.angle - panel.angle for panel in panels ]) 

        # Solve system.
        vts = At*strengths + bt

        # Update panel velocities and pressure coefficients.
        for (panel, vt) in zip(panels, vts)
            panel.vt = vt
            panel.cp = pressure ? pressureCoefficient2D(0., vt, uniform.magnitude) : nothing
        end

        return vts
    end

    function liftCoefficient(gamma :: Float64, panels :: Array{VortexSourcePanel2D, 1}, uniform :: Uniform2D)
        # Lift coefficient
        xs = [ panel.xs for panel in panels ]
        c = abs(maximum(xs) - minimum(xs))
        cl = gamma/(0.5*uniform.magnitude*c)*sum([ panel.length for panel in panels ])

        return cl
    end

    # Compute errors
    errorCalc(panels :: Array{<: Panel, 1}) = sum([ panel.strength*panel.length for panel in panels ])

    # Compute velocities
    velocity(panels :: Array{<: Panel, 1}, uniform :: Uniform2D, x :: Float64, y :: Float64) = uniform.velocity(x,y) .+ foldl( (u, v) -> u .+ v, [ panel.velocity(panel.strength, x, y) for panel in panels ] )
    
    # Compute potential
    potential(panels :: Array{<: Panel, 1}, uniform :: Uniform2D, x :: Float64, y :: Float64) = uniform.potential(x, y) + sum([ panel.potential(x, y) for panel in panels ])

    # gridCompute(xs :: Array{<: Real}) = [ velocity(..., x...) for x in xs ]

    # Euclidean norm of a vector
    mag(x) = sqrt(sum(s -> s^2, x))

    # Infinite vortices
    infiniteVortices(strength :: Float64, spacing :: Float64, x :: Float64, y :: Float64) = (strength/(2*spacing)*sinh(2π*y/spacing)/(cosh(2π*y/spacing) - cos(2π*x/spacing)), strength/(2*spacing)*sin(2π*x/spacing)/(cosh(2π*y/spacing) - cos(2π*x/spacing)))

    # Joukowski transformation
    joukowski(z :: Complex, c :: Float64) = z + c^2/z

    # Pressure coefficient
    pressureCoefficient2D(u :: Float64, v :: Float64, freestream_speed :: Float64) = 1. - (u^2 + v^2)/freestream_speed^2

    function cosinePanels(x :: Array{<: Real}, y :: Array{<: Real}, n :: Integer = 40)
        """
        Discretises a geometry consisting of x and y coordinates into panels by projecting the x-coordinate of a circle onto the geometry.
        """
        r = (maximum(x) - minimum(x))/2
        x_center = (maximum(x) + minimum(x))/2
        x_ends = x_center .+ r*cos.(range(0.0, stop = 2π, length = n+1))
        y_ends = zeros(length(x_ends))
    
        x, y = push!(x, x[1]), push!(y, y[1])
    
        j = 1
        for i in 1:n
            while j < length(x)
                if (x[j] <= x_ends[i] <= x[j+1]) | (x[j+1] <= x_ends[i] <= x[j])
                    break
                else
                    j += 1
                end
            end
            m = (y[j+1] - y[j])/(x[j+1] - x[j])
            c = y[j+1] - m*x[j+1]
            y_ends[i] = m*x_ends[i] + c
        y_ends[n + 1] = y_ends[1]
        end

        panels = [ Panel2D(x_ends[i], y_ends[i], x_ends[i+1], y_ends[i+1]) for i in 1:n ]

        return panels
    end
end
