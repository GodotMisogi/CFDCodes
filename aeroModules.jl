module panelSolver

using QuadGK

abstract type Solution end

struct Source2D <: Solution
    strength
    x_0
    y_0
    velocity
    stream
    potential
    function Source2D(strength, x_0, y_0)
        velocity(x, y) = (strength/(2π)*(x - x_0)/((x - x_0)^2 + (y - y_0)^2), strength/(2π)*(y - y_0)/((x - x_0)^2 + (y - y_0)^2))
        stream(x, y) = strength/(2π)*atan(y - y_0, x - x_0)
        potential(x, y) = strength/(4π)*log((x - x_0)^2 + (y - y_0)^2)
        new(strength, x_0, y_0, velocity, stream, potential)
    end
end

struct Uniform2D <: Solution
    magnitude
    angle
    deg
    velocity
    stream
    potential
    function Uniform2D(magnitude, deg)
        angle = deg2rad(deg)
        
        velocity(x, y) = (magnitude*cos(angle), magnitude*sin(angle))
        stream(x, y) = magnitude*(y*cos(angle) - x*sin(angle))
        potential(x, y) = magnitude*(x*cos(angle) + y*sin(angle))
        new(magnitude, angle, deg, velocity, stream, potential)
    end
end

struct Doublet2D <: Solution
    strength
    x_0
    y_0
    velocity
    stream
    potential
    function Doublet2D(strength, x_0, y_0)
        velocity(x, y) = (-strength/(2π)*((x - x_0)^2 - (y - y_0)^2)/((x - x_0)^2 + (y - y_0)^2)^2, -strength/(2π)*2*(x - x_0)*(y - y_0)/((x - x_0)^2 + (y - y_0)^2)^2)
        stream(x, y) = -strength/(2π)*(y - y_0)/((x - x_0)^2 + (y - y_0)^2)
        potential(x, y) = -strength/(2π)*(x - x_0)/((x - x_0)^2 + (y - y_0)^2)
        new(strength, x_0, y_0, velocity, stream, potential)
    end
end

struct Vortex2D <: Solution
    strength
    x_0
    y_0
    velocity
    stream
    potential
    function Vortex2D(strength, x_0, y_0)
        velocity(x, y) = (-strength/(2π)*(y - y_0)/((x - x_0)^2 + (y - y_0)^2), strength/(2π)*(x - x_0)/((x - x_0)^2 + (y - y_0)^2))
        stream(x, y) = -strength/(4π)*log((x - x_0)^2 + (y - y_0)^2)
        potential(x, y) = strength/(2π)*atan(y - y_0, x - x_0)
        new(strength, x_0, y_0, velocity, stream, potential)
    end
end

mutable struct SourcePanel2D <: Solution
    xs
    ys
    xe 
    ye
    xc
    yc
    length
    angle
    loc :: T where T <: AbstractString
    strength
    vt
    cp
    integral
    potential
    velocity
    function SourcePanel2D(xs, ys, xe, ye, strength=0.0, vt=0.0, cp=0.0)
        xc, yc = (xe + xs)/2.0, (ye + ys)/2.0
        length = mag([ xe - xs, ye - ys ])
        angle = xe <= xs ? acos((ye - ys)/length) : π + acos(-(ye - ys)/length)
        loc = angle <= π ? "upper" : "lower"

        integral(x, y, dxdz, dydz) = quadgk(s -> ((x - (xs - s*sin(angle)))*dxdz + (y - (ys + s*cos(angle)))*dydz)/((x - (xs - s*sin(angle)))^2 + (y - (ys + s*cos(angle)))^2), 0., length, atol=1e-10)[1]
        potential(x, y) = quadgk( (s) -> log((x - (xs - s*sin(angle)))^2 + (y - (ys + s*cos(angle)))^2), 0, length, atol=1e-10)[1]
        velocity(strength, x, y) = strength/(2π).*(integral(x, y, 1., 0.), integral(x, y, 0., 1.))
        new(xs, ys, xe, ye, xc, yc, length, angle, loc, strength, vt, cp, integral, potential, velocity)
    end 
end

struct SourcePanelSolver2D <: Solution
    """
    Applies the source panel method to a list of panels and a uniform flow.
    """
    panels :: Array{SourcePanel2D, 1}
    uniform :: Uniform2D
    error
    velocity
    potential
    solveStrengths
    aerodynamicsss
    function SourcePanelSolver2D(panels, uniform)
        function solveStrengths()
            """
            Solves for the source strengths of all panels.
            """
            # Construct matrix for normal direction.
            An = [ i == j ? 0.5 : 0.5/π*panel_i.integral(panel_j.xc, panel_j.yc, cos(panel_j.angle), sin(panel_j.angle)) for (j, panel_j) in enumerate(panels), (i, panel_i) in enumerate(panels) ]

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
        
        function aerodynamicsss(pressure = true)
            """
            Solves for the velocities of all panels and their pressure coefficients.
            """
            # Construct matrix for tangential direction.
            At = [ i == j ? 0. : 0.5/π*panel_i.integral(panel_j.xc, panel_j.yc, -sin(panel_j.angle), cos(panel_j.angle)) for (j, panel_j) in enumerate(panels), (i, panel_i) in enumerate(panels) ]

            # Construct freestream RHS.
            bt = @. uniform.magnitude*sin([ uniform.angle - panel.angle for panel in panels ]) 

            # Solve system.
            vts = At*[ panel.strength for panel in panels ] + bt
    
            # Update panel velocities and pressure coefficients.
            for (panel, vt) in zip(panels, vts)
                panel.vt = vt
                panel.cp = pressure ? pressureCoefficient2D(0., vt, uniform.magnitude) : 0
            end
    
            return vts
        end
        # Compute errors
        error() = sum([ panel.strength*panel.length for panel in panels ])
        # Compute velocities
        velocity(x, y) = uniform.velocity(x,y) .+ foldl( (u, v) -> u .+ v, [ panel.velocity(panel.strength, x, y) for panel in panels ] )
        # Compute potential
        potential(x, y) = uniform.potential(x, y) .+ sum([ panel.potential(x, y) for panel in panels ])
        new(panels, uniform, error, velocity, potential, solveStrengths, aerodynamicsss)
    end
end

function panelSolver2D(aeroproblem :: SourcePanelSolver2D, uniform :: Uniform2D)
    strengths = aeroproblem.solveStrengths()
    vts = aeroproblem.aerodynamicsss()
    error = aeroproblem.error()
    return strengths, vts, error
end


mutable struct VortexPanel2D <: Solution
    xs
    ys
    xe 
    ye
    x_14
    y_14
    xc
    yc
    length
    angle
    loc :: T where T <: AbstractString
    strength
    vt
    cp
    influence
    # integral
    potential
    velocity
    # function VortexPanel2D(xs, ys, xe, ye, strength=0.0, vt=0.0, cp=0.0)
    #     xc, yc = (xe + xs)/2.0, (ye + ys)/2.0
    #     length = mag([ xe - xs, ye - ys ])
    #     angle = xe <= xs ? atan(ye - ys, xe - xs) : π + atan(-ye - ys, xe - xs)
    #     loc = angle <= π ? "upper" : "lower"

    #     integral(x, y, dxdz, dydz) = quadgk(s -> ((x - (xs - s*sin(angle)))*dydz - (y - (ys + s*cos(angle)))*dxdz)/((x - (xs - s*sin(angle)))^2 + (y - (ys + s*cos(angle)))^2), 0., length, atol=1e-10)[1]
    #     potential(x, y) = quadgk( (s) -> atan((y - (ys + s*cos(angle)))/(x - (xs - s*sin(angle)))), 0, length, atol=1e-10)[1]
    #     velocity(strength, x, y) = strength/(2π).*(integral(x, y, 0., 1.), integral(x, y, -1., 0.))
    #     new(xs, ys, xe, ye, xc, yc, length, angle, loc, strength, vt, cp, integral, potential, velocity)
    # end 
    function VortexPanel2D(xs, ys, xe, ye, strength=0.0, vt=0.0, cp=0.0)
        x_14, y_14 = xs + (xe - xs)/4., ys + (ye - ys)/4.
        xc, yc = xs + (3/4)*(xe - xs), ys + (3/4)*(ye - ys)
        # xc, yc = (xe + xs)/2.0, (ye + ys)/2.0
        length = mag([ xe - xs, ye - ys ])
        angle = atan(ye - ys, xe - xs)
        loc = angle <= π ? "upper" : "lower"

        influence(x, y, dxdz, dydz) = 1/(2π*((x - x_14)^2 + (y - y_14)^2))*((y - y_14)*dxdz - (x - x_14)*dydz)
        potential(x, y) = strength/(2π)*atan(y - y_14, x - x_14)
        velocity(x, y) = strength/(2π*((x - x_14)^2 + (y - y_14)^2)).*(y - y_14, -(x - x_14))
        new(xs, ys, xe, ye, x_14, y_14, xc, yc, length, angle, loc, strength, vt, cp, influence, potential, velocity)
    end 
end

struct VortexPanelSolver2D <: Solution
    panels :: Array{VortexPanel2D, 1}
    uniform :: Uniform2D
    error
    velocity
    potential
    solveStrengths
    aerodynamicsss
    liftCoefficient
    function VortexPanelSolver2D(panels, uniform)
        function solveStrengths()
            """
            Solves for the source strengths of all panels.
            """
            # Construct matrix for normal direction.
            An = [ panel_i.influence(panel_j.xc, panel_j.yc, sin(panel_j.angle), cos(panel_j.angle)) for (j, panel_j) in enumerate(panels), (i, panel_i) in enumerate(panels) ]

            # Construct freestream RHS.
            bn = @. -uniform.magnitude*sin([ uniform.angle + panel.angle for panel in panels ]) 
            # Solve system.
            sources = An\bn
    
            # Update panel strengths.
            for (panel, strength) in zip(panels, sources)
                panel.strength = strength
            end

            return sources
        end
        
        function aerodynamicsss(pressure = true)
            """
            Solves for the velocities of all panels and their pressure coefficients.
            """
            # Construct matrix for tangential direction.
            # At = [ panel_i.influence(panel_j.x_34, panel_j.y_34, -sin(panel_j.angle), cos(panel_j.angle)) for (j, panel_j) in enumerate(panels), (i, panel_i) in enumerate(panels) ]

            # Construct freestream RHS.
            # bt = @. uniform.magnitude*cos([ uniform.angle + panel.angle for panel in panels ]) 
    
            # Solve system.
            # vts = At*[ panel.strength for panel in panels ] + bt
    
            # Update panel velocities and pressure coefficients.
            for panel in panels
                # panel.vt = 0
                panel.cp = pressure ? panel.strength/(uniform.magnitude*panel.length) : 0
            end
    
            return cp
        end
        function liftCoefficient()
            # Lift coefficient
            xs = [ panel.xs for panel in panels ]
            c = abs(maximum(xs) - minimum(xs))
            cl = 1.0/(0.5*uniform.magnitude*c)*sum([ panel.strength for panel in panels ])

            return cl
        end
        # Compute errors
        error() = sum([ panel.strength*panel.length for panel in panels ])
        # Compute velocities
        velocity(x, y) = uniform.velocity(x,y) .+ foldl( (u, v) -> u .+ v, [ panel.velocity(x, y) for panel in panels ] )
        zip
        # Compute potential
        potential(x, y) = uniform.potential(x, y) .+ sum([ panel.potential(x, y) for panel in panels ])
        new(panels, uniform, error, velocity, potential, solveStrengths, aerodynamicsss, liftCoefficient)
    end
end

function panelSolver2D(aeroproblem :: VortexPanelSolver2D, uniform :: Uniform2D)
    strengths = aeroproblem.solveStrengths()
    vts = aeroproblem.aerodynamicsss()
    cl = aeroproblem.liftCoefficient()
    error = aeroproblem.error()
    return strengths, vts, cl, error
end

mutable struct DoubletPanel2D <: Solution
    xs
    ys
    xe 
    ye
    xc
    yc
    length
    angle
    loc :: T where T <: AbstractString
    normal
    tangent
    strength
    vt
    cp
    potential
    velocity
    function DoubletPanel2D(xs, ys, xe, ye, strength=0.0, vt=0.0, cp=0.0)
        length = mag([ xe - xs, ye - ys ])
        angle = atan(ye - ys, xe - xs)
        loc = (π/2 <= angle <= π) || (-π <= angle <=-π/2) ? "lower" : "upper"
        normal = (-sin(angle), cos(angle))
        tangent = (cos(angle), sin(angle))
        
        # Collocation points in global coordinates
        xc, yc = (xe + xs)/2.0, (ye + ys)/2.0

        # Panel coordinates
        xsl, ysl, xel, yel = 0., 0., length, 0.
        xcl, ycl = length/2., 0

        potential(strength, x, y) = strength/(4π)*(x - xc)/((x - xc)^2 + (y - yc)^2)
        function velocity(strength, xg, yg) 
            x, y = panelCoords(xg, yg, xs, ys, angle)
            return invRotation(strength/(2π).*((y/((x - xsl)^2 + y^2) - y/((x - xel)^2 + y^2)), -1*((x - xsl)/((x - xsl)^2 + y^2) - (x - xel)/((x - xel)^2 + y^2)))..., angle)
        end
        new(xs, ys, xe, ye, xc, yc, length, angle, loc, normal, tangent, strength, vt, cp, potential, velocity)
    end 
end

struct DoubletPanelSolver2D <: Solution
    """
    Applies the doublet panel method to a list of panels and a uniform flow.
    """
    panels :: Array{DoubletPanel2D, 1}
    uniform :: Uniform2D
    error
    velocity
    potential
    solveStrengths
    aerodynamicsss
    function DoubletPanelSolver2D(panels, uniform)
        num_panels = length(panels)
        woke_panel = DoubletPanel2D(panels[end].xe, panels[end].ye, 1000*panels[end].xe, panels[end].ye)
        u = uniform.magnitude.*(cos(uniform.angle), sin(uniform.angle))
        function solveStrengths()
            """
            Solves for the doublet strengths of all panels.
            """

            # Construct doublet influence matrix for normal direction.
            doublet = [ i == j ? -2/(π*panel_i.length) : sum(panel_j.velocity(1., panel_i.xc, panel_i.yc).*panel_i.normal) for (i, panel_i) in enumerate(panels), (j, panel_j) in enumerate(panels) ]

            woke = [ sum(woke_panel.velocity(1., panel.xc, panel.yc).*panel.normal) for panel in panels ]

            # Kutta condition
            kutta = zeros(num_panels + 1)
            # kutta[1] = 1.
            # kutta[end-1] = -1.
            # kutta[end] = 1.

            # Stronger Kutta condition
            kutta[1] = 1.
            kutta[2] = -1.
            kutta[end-2] = 1.
            kutta[end-1] = -1.

            # Matrix construction
            An = zeros(num_panels + 1, num_panels + 1)
            An[1:end-1, 1:end-1] = doublet
            An[1:end-1, end] = woke
            An[end, :] = kutta
            # print(An)

            # Construct freestream RHS.
            bn = push!([ -sum(u.*(panel.normal)) for panel in panels], 0)

            # Solve system.
            strengths = An\bn

            # Update panel strengths.
            for (panel, strength) in zip(panels, strengths)
                panel.strength = strength
            end

            return strengths
        end
        
        function aerodynamicsss(pressure = true)
            """
            Solves for the velocities of all panels and their pressure coefficients.
            """
            push!(panels, woke_panel)
            
            # Update panel velocities, pressure and lift coefficients.
            cl = 0
            for (i, panel) in enumerate(panels)
                if panel == panels[1]
                    R = mag([panels[2].xc - panels[1].xc, panels[2].yc - panels[1].yc])
                    panel.vt = (panels[2].strength - panels[1].strength)/R + sum(u.*(panel.tangent))
                elseif panel == panels[end]
                    R = mag([panels[end].xc - panels[end-1].xc, panels[end].yc - panels[end-1].yc])
                    panel.vt = (panels[end].strength - panels[end-1].strength)/R + sum(u.*(panel.tangent))
                else
                    R = mag([panels[i + 1].xc - panels[i - 1].xc, panels[i + 1].yc - panels[i - 1].yc])
                    panel.vt = (panels[i + 1].strength - panels[i - 1].strength)/R + sum(u.*(panel.tangent))
                end
                panel.cp = pressure ? pressureCoefficient2D(0., panel.vt, uniform.magnitude) : 0.
                cl = cl - panel.vt*R
            end

            cl = cl/uniform.magnitude

            panels = push!(panels, woke_panel)
            xs = [ panel.xs for panel in panels ]
            c = abs(maximum(xs) - minimum(xs))
            vts = [ panel.vt for panel in panels ]

            return (vts, cl)
        end
        
        # Compute errors
        error() = sum([ panel.strength*panel.length for panel in panels ])
        # Compute velocities
        velocity(x, y) = uniform.velocity(x, y) .+ foldl( (u, v) -> u .+ v, [ panel.velocity(panel.strength, x, y) for panel in panels ] )
        # Compute potential 
        potential(x, y) = uniform.potential(x, y) .+ sum([ panel.potential(panel.strength, x, y) for panel in panels ])
        new(panels, uniform, error, velocity, potential, solveStrengths, aerodynamicsss)
    end
end

function panelSolver2D(aeroproblem :: DoubletPanelSolver2D, uniform :: Uniform2D)
    strengths = aeroproblem.solveStrengths()
    vts, cl = aeroproblem.aerodynamicsss()
    error = aeroproblem.error()
    return strengths, vts, cl, error
end

struct VortexSourcePanelSolver2D <: Solution
    """
    Applies the vortex source panel method to a list of panels and a uniform flow.
    """
    panels :: Array{SourcePanel2D, 1}
    uniform :: Uniform2D
    gamma
    error
    velocity
    potential
    solveStrengths
    aerodynamicsss
    liftCoefficient
    function VortexSourcePanelSolver2D(panels, uniform, gamma=0.0)            
        # Construct source matrix for normal direction.
        source_matrix = [ i == j ? 0.5 : 0.5/π*panel_i.integral(panel_j.xc, panel_j.yc, cos(panel_j.angle), sin(panel_j.angle)) for (j, panel_j) in enumerate(panels), (i, panel_i) in enumerate(panels) ]

        # Vortex matrix
        vortex_matrix = [ i == j ? 0.0 : -0.5/π*panel_i.integral(panel_j.xc, panel_j.yc, sin(panel_j.angle), -cos(panel_j.angle)) for (j, panel_j) in enumerate(panels), (i, panel_i) in enumerate(panels) ]

        num_panels = length(panels)

        function solveStrengths()
            """
            Solves for the source strengths of all panels.
            """
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
            bn = push!( -uniform.magnitude .* cos.([ uniform.angle - panel.angle for panel in panels ]), -uniform.magnitude * (sin(uniform.angle - panels[1].angle) + sin(uniform.angle - panels[end].angle)))

            # Solve system.
            strengths = An\bn

            # Update panel strengths.
            for (panel, strength) in zip(panels, strengths)
                panel.strength = strength
            end
            
            # Update vorticity.
            gamma = strengths[end]

            return strengths
        end
        function aerodynamicsss(pressure=true)
            """
            Solves for the velocities of all panels and their pressure coefficients.
            """
            At = zeros(num_panels, num_panels + 1)

            # Construct matrix for tangential direction.
            At[:, 1:end-1] = vortex_matrix
            At[:, end] = -sum(source_matrix, dims=2)
            
            # Construct freestream RHS.
            bt = @. uniform.magnitude*sin([ uniform.angle - panel.angle for panel in panels ]) 

            # Solve system.
            vts = At*push!([ panel.strength for panel in panels ], gamma) + bt

            # Update panel velocities and pressure coefficients.
            for (panel, vt) in zip(panels, vts)
                panel.vt = vt
                panel.cp = pressure ? pressureCoefficient2D(0., vt, uniform.magnitude) : 0
            end

            return vts
        end
        function liftCoefficient()
            # Lift coefficient
            xs = [ panel.xs for panel in panels ]
            c = abs(maximum(xs) - minimum(xs))
            cl = gamma/(0.5*uniform.magnitude*c)*sum([ panel.length for panel in panels ])

            return cl
        end
        # Compute errors
        error() = sum([ panel.strength*panel.length for panel in panels ])
        # Compute velocities
        velocity(x, y) = uniform.velocity(x,y) .+ foldl( (u, v) -> u .+ v, [ panel.velocity(panel.strength, x, y) for panel in panels ] )
        # Compute potential
        potential(x, y) = uniform.potential(x, y) .+ sum([ panel.potential(x, y) for panel in panels ])
        new(panels, uniform, gamma, error, velocity, potential, solveStrengths, aerodynamicsss, liftCoefficient)
    end
end

function panelSolver2D(aeroproblem :: VortexSourcePanelSolver2D, uniform :: Uniform2D)
    strengths = aeroproblem.solveStrengths()
    vts = aeroproblem.aerodynamicsss()
    cl = aeroproblem.liftCoefficient()
    error = aeroproblem.error()
    return strengths, vts, cl, error
end

mutable struct DoubletSourcePanel2D <: Solution
    xs
    ys
    xe 
    ye
    xc
    yc
    length
    angle
    loc :: T where T <: AbstractString
    normal
    tangent
    doublet_strength
    source_strength
    vt
    cp
    sourceInfluence
    doubletInfluence
    potential
    velocity
    function DoubletSourcePanel2D(xs, ys, xe, ye, doublet_strength=0.0, source_strength=1.0, vt=0.0, cp=0.0)
        # Properties
        length = mag([xe - xs, ye - ys])
        angle = atan(ye - ys, xe - xs)
        loc = (π/2 <= angle <= π) || (-π <= angle <=-π/2) ? "lower" : "upper"
        normal = (-sin(angle), cos(angle))
        tangent = (cos(angle), sin(angle))

        # Collocation points
        xc, yc = (xe + xs)/2.0, (ye + ys)/2.0

        # Panel coordinates
        xsl, ysl, xel, yel = 0., 0., length, 0.
        xcl, ycl = length/2., 0

        function sourceInfluence(xg, yg) 
            x, y = panelCoords(xg, yg, xs, ys, angle)
            return 1/(4π)*((x - xsl)*log((x - xsl)^2 + y^2) - (x - xel)log((x - xel)^2 + y^2) + 2y*(atan(y, x - xel) - atan(y, x)))
        end
        function doubletInfluence(xg, yg)
            x, y = panelCoords(xg, yg, xs, ys, angle) 
            return -1/(2π)*(atan(y, x - xel) - atan(y, x - xsl))
        end
        
        function potential(doublet_strength, source_strength, xg, yg)
            return doublet_strength.*doubletInfluence(xg, yg) .+ source_strength.*sourceInfluence(xg, yg)
        end
        
        function velocity(doublet_strength, source_strength, xg, yg) 
            x, y = panelCoords(xg, yg, xs, ys, angle)
            return invRotation(doublet_strength/(2π).*((y/((x - xsl)^2 + y^2) - y/((x - xel)^2 + y^2)), -((x - xsl)/((x - xsl)^2 + y^2) - (x - xel)/((x - xel)^2 + y^2))) .+ source_strength/(4π).*(log(((x - xsl)^2 + y^2)/((x - xel)^2 + y^2)), 2*(atan(y, x - xel) - atan(y, x - xsl)))..., angle)
        end        
        new(xs, ys, xe, ye, xc, yc, length, angle, loc, normal, tangent, doublet_strength, source_strength, vt, cp, sourceInfluence, doubletInfluence, potential, velocity)
    end 
end

struct DoubletSourcePanelSolver2D <: Solution
    """
    Applies the constant-potential doublet-source panel method to a list of panels and a uniform flow.
    """
    panels :: Array{DoubletSourcePanel2D, 1}
    uniform :: Uniform2D
    error
    velocity
    potential
    solveStrengths
    aerodynamicsss
    liftCoefficient
    function DoubletSourcePanelSolver2D(panels, uniform, sources=false)
        num_panels = length(panels)
        woke_panel = DoubletSourcePanel2D(panels[end].xe, panels[end].ye, 100000*panels[end].xe, panels[end].ye)
        u = uniform.magnitude.*(cos(uniform.angle), sin(uniform.angle))
        function solveStrengths()
            """
            Solves for the strengths of all panels.
            """
            # Wake vector
            woke_vector = [ woke_panel.doubletInfluence(panel.xc, panel.yc) for panel in panels ]

            # Doublet matrix
            doublet_matrix = [ i == j ? 0.5 : panel_j.doubletInfluence(panel_i.xc, panel_i.yc) for (i, panel_i) in enumerate(panels), (j, panel_j) in enumerate(panels) ]
            
            # Kutta condition
            # kutta = zeros(num_panels + 1)
            # kutta[1] = 1.
            # kutta[end-1] = -1.
            # kutta[end] = 1.

            # Stronger Kutta condition
            # kutta = zeros(num_panels + 1)
            # kutta[1] = 1.
            # kutta[2] = -1.
            # kutta[end-2] = 1.
            # kutta[end-1] = -1.

            # LHS with explicit Kutta condition
            # An = zeros(num_panels + 1, num_panels + 1)
            # An[1:end-1, 1:end-1] = doublet_matrix
            # An[1:end-1, end] = woke_vector 
            # An[end, :] = kutta

            # LHS with implicit Kutta condition
            An = zeros(num_panels, num_panels)
            An = doublet_matrix
            An[:, 1] = doublet_matrix[:, 1] - woke_vector 
            An[:, end] = doublet_matrix[:, end] + woke_vector

            b = zeros(num_panels)
            if sources
                # Source influence matrix
                source_matrix = [ panel_j.sourceInfluence(panel_i.xc, panel_i.yc) for (i, panel_i) in enumerate(panels), (j, panel_j) in enumerate(panels) ]

                # Source vector (σ = ∂Φ/∂n · n)
                source_vector = [ sum(u.*(panel.normal)) for panel in panels ]

                # Update panel source strengths
                for (panel, source) in zip(panels, source_vector)
                    panel.source_strength = source
                end

                # Source RHS
                b[1:end] = -source_matrix*source_vector
            else
                # Freestream RHS
                b[1:end] = [ -sum(u.*(panel.xc, panel.yc)) for panel in panels ]
            end

            # Solve system
            strengths = An\b

            # Update panel doublet strengths
            for (panel, strength) in zip(panels, strengths)
                panel.doublet_strength = strength
            end

            woke_panel.doublet_strength = strengths[end] - strengths[1]

            return strengths
        end
        function aerodynamicsss(pressure=true)
            """
            Solves for the velocities of all panels and their pressure coefficients.
            """
            
            xs = [ panel.xs for panel in panels ]
            c = abs(maximum(xs) - minimum(xs))

            # Update panel velocities, pressure and lift coefficients.
            cl = 0
            if sources 
                phi = [ sum(u.*(panel.xc, panel.yc)) + panel.doublet_strength for panel in panels ]
                colengths = [ mag([panel1.xc - panel2.xc, panel1.yc - panel2.yc]) for (panel1, panel2) in zip(panels[1:end-1], panels[2:end]) ]
                vts = [ (phi1 - phi2)/r for (phi1, phi2, r) in zip(phi[1:end-1], phi[2:end], colengths) ]
                for (vt, panel) in zip(vts, panels)
                    panel.vt = vt
                    panel.cp = pressure ? pressureCoefficient2D(0., vt, uniform.magnitude) : 0.
                end
                cl = (sum([ -panel.cp*r*cos(panel.angle)/c for (panel, r) in zip(panels[1:end-1], colengths) ]), 2*woke_panel.doublet_strength/uniform.magnitude)
            else
                for (i, panel) in enumerate(panels)
                    if panel == panels[1]
                        R = mag([panels[2].xc - panels[1].xc, panels[2].yc - panels[1].yc])
                        panel.vt = (panels[2].doublet_strength - panels[1].doublet_strength)/R
                    elseif panel == panels[end]
                        R = mag([panels[end].xc - panels[end-1].xc, panels[end].yc - panels[end-1].yc])
                        panel.vt = (panels[end].doublet_strength - panels[end-1].doublet_strength)/R
                    else
                        R = mag([panels[i+1].xc - panels[i-1].xc, panels[i+1].yc - panels[i-1].yc])
                        panel.vt = (panels[i+1].doublet_strength - panels[i-1].doublet_strength)/R
                    end
                    panel.cp = pressure ? pressureCoefficient2D(0., panel.vt, uniform.magnitude) : 0.
                    cl = cl - panel.vt*R
                end
                cl = (cl/uniform.magnitude, -2*woke_panel.doublet_strength/uniform.magnitude)
                vts = [ panel.vt for panel in panels ]
            end

            return (vts, cl)
        end
        # Compute errors
        error() = sum([ panel.doublet_strength*panel.length for panel in panels])
        # Compute velocities
        velocity(x, y) = uniform.velocity(x, y) .+ foldl( (u, v) -> u .+ v, [ panel.velocity(panel.doublet_strength, panel.source_strength, x, y) for panel in panels ] )
        # Compute potential
        potential(x, y) = uniform.potential(x, y) .+ sum([ panel.potential(panel.doublet_strength, panel.source_strength, x, y) for panel in panels ])
        new(panels, uniform, error, velocity, potential, solveStrengths, aerodynamicsss)
    end
end

function panelSolver2D(aeroproblem :: DoubletSourcePanelSolver2D, uniform :: Uniform2D)
    strengths = aeroproblem.solveStrengths()
    vts, cl = aeroproblem.aerodynamicsss()
    error = aeroproblem.error()
    return strengths, vts, cl, error
end

# Performs velocity and potential calculations on a grid
function gridData(objects :: Array{<:Solution, 1}, xs, pressure = true)
    vels = foldl( (v1, v2) -> [ u .+ v for (u, v) in zip(v1, v2)], [ velocity(object, xs) for object in objects ])
    pots = foldl( (v1, v2) -> v1 .+ v2, [ potential(object, xs) for object in objects ])
end

function gridData(object :: Solution, xs)
    vels = velocity(object, xs)
    pots = potential(object, xs)
    
    return vels, pots
end

# Performs velocity and potential computations for an object on a grid
velocity(object :: Solution, xs) = [ object.velocity(x...) for x in xs ] 
potential(object :: Solution, xs) = [ object.potential(x...) for x in xs ]

# Transforms (x, y) to the coordinate system with (x_s, y_s) as origin oriented at angle_s.
panelCoords(x, y, x_s, y_s, angle_s) = rotation(x - x_s, y - y_s, angle_s)

# Rotation matrices
invRotation(x, y, angle) = (x*cos(angle) - y*sin(angle), x*sin(angle) + y*cos(angle))
rotation(x, y, angle) = (x*cos(angle) + y*sin(angle), -x*sin(angle) + y*cos(angle))

# Pressure coefficient
pressureCoefficient2D(vels :: Array, freestream) = pressureCoefficient2D.(first.(vels), last.(vels), freestream)
pressureCoefficient2D(u, v, freestream_speed) = 1. - mag([u,v])^2/freestream_speed^2

# Euclidean norm of a vector
mag(x) = sqrt(sum(s -> s^2, x))

# Infinite vortices
infiniteVortices(strength :: Float64, spacing :: Float64, x :: Float64, y :: Float64) = (strength/(2*spacing)*sinh(2π*y/spacing)/(cosh(2π*y/spacing) - cos(2π*x/spacing)), strength/(2*spacing)*sin(2π*x/spacing)/(cosh(2π*y/spacing) - cos(2π*x/spacing)))

# Joukowski transformation
joukowski(z :: Complex, c) = z + c^2/z

function cosinePanels(x :: Array{<:Real}, y :: Array{<:Real}, n = 40)
    """
    Discretises a geometry consisting of x and y coordinates into panels by projecting the x-coordinate of a circle onto the geometry.
    """
    r = (maximum(x) - minimum(x))/2.
    x_center = (maximum(x) + minimum(x))/2.
    if any(y -> y < 0., y)
        x_circ = x_center .+ r.*cos.(range(-π, stop = 0.0, length = n + 1))
    else
        x_circ = x_center .+ r.*cos.(range(0.0, stop = π, length = n + 1))
    end

    x_ends = copy(x_circ)
    y_ends = zeros(length(x_ends))

    x, y = push!(x, x[1]), push!(y, y[1])

    j = 1
    for i in 1:n+1
        while j < length(x)
            if ((x[j] <= x_ends[i] <= x[j+1]) | (x[j+1] <= x_ends[i] <= x[j]))
                break
            else
                j += 1
            end
        end
        m = (y[j+1] - y[j])/(x[j+1] - x[j])
        c = y[j+1] - m*x[j+1]
        y_ends[i] = m*x_ends[i] + c
        # y_ends[n+1] = y_ends[n]
    end

    return (x_ends, y_ends)
end


function cosineAirfoil(x :: Array{<:Real}, y :: Array{<:Real}, n :: Integer = 40)
    """
    Discretises a geometry consisting of x and y coordinates into panels by projecting the x-coordinate of a circle onto the geometry.
    """
    r = (maximum(x) - minimum(x))/2.
    x_center = (maximum(x) + minimum(x))/2.
    x_circ = x_center .+ r.*cos.(range(0.0, stop = 2π, length = n+1))

    x_ends = copy(x_circ)
    y_ends = zeros(length(x_ends))

    x, y = push!(x, x[1]), push!(y, y[1])

    j = 1
    for i in 1:n
        while j < length(x)
            if ((x[j] <= x_ends[i] <= x[j+1]) | (x[j+1] <= x_ends[i] <= x[j]))
                break
            else
                j += 1
            end
        end
        m = (y[j+1] - y[j])/(x[j+1] - x[j])
        c = y[j+1] - m*x[j+1]
        y_ends[i] = m*x_ends[i] + c
        y_ends[n+1] = y_ends[1]
    end

    return (x_ends, y_ends)
end

function NACA4(digits :: Tuple, chord, n, closed_te=false)
    
    # Airfoil characteristics
    camber = digits[1]/100
    position = digits[2]/10
    t_by_c = (10*digits[3] + digits[4])/100

    # Cosine spacing
    angles = range(0, stop = π, length = n + 1)
    xs = reverse([ chord*(1 - 0.5*(1 - cos(beta))) for beta in angles ] )

    # Thickness distribution
    thickness = closed_te ? [ 5*t_by_c*chord*(0.2969*sqrt(xc/chord) - 0.126*xc/chord - 0.3516*(xc/chord)^2 + 0.2843*(xc/chord)^3 - 0.1036*(xc/chord)^4) for xc in xs ] : [ 5*t_by_c*chord*(0.2969*sqrt(xc/chord) - 0.126*xc/chord - 0.3516*(xc/chord)^2 + 0.2843*(xc/chord)^3 - 0.1015*(xc/chord)^4) for xc in xs ] 
    
    if position == 0.0 || camber == 0.0
        x_upper = xs
        y_upper = thickness
        x_lower = xs
        y_lower = -thickness
    else
        camberline = [ xc <= position*camber ? (camber/position^2)*xc*(2*position-xc/chord) : (camber/(1-position)^2)*(chord - xc)*(1 + xc/chord - 2*position) for xc in xs ]
        gradients = [ xc <= position*camber ? atan((2*camber/position)*(-xc/(chord*position) + 1)) : atan((2*camber/(1 - position^2))*(position - xc/chord)) for xc in xs ] 
        x_upper = [ xc - thicc*sin(thot) for (xc, thicc, thot) in zip(xs, thickness, gradients) ] 
        y_upper = [ cam + thicc*cos(thot) for (cam, thicc, thot) in zip(camberline, thickness, gradients) ]
        x_lower = [ xc + thicc*sin(thot) for (xc, thicc, thot) in zip(xs, thickness, gradients) ] 
        y_lower = [ cam - thicc*cos(thot) for (cam, thicc, thot) in zip(camberline, thickness, gradients) ]
    end

    (X, Y) = append!(reverse(x_upper), x_lower), append!(reverse(y_upper), y_lower)
        
    return (X, Y)
end

end