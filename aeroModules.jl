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
        tangentialVelocities
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
            
            function tangentialVelocities(pressure = true)
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
            new(panels, uniform, error, velocity, potential, solveStrengths, tangentialVelocities)
        end
    end

    function panelSolver2D(aeroproblem :: SourcePanelSolver2D, uniform :: Uniform2D)
        strengths = aeroproblem.solveStrengths()
        vts = aeroproblem.tangentialVelocities()
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
        x_34
        y_34
        # xc
        # yc
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
            x_34, y_34 = xs + (3/4)*(xe - xs), ys + (3/4)*(ye - ys)
            # xc, yc = (xe + xs)/2.0, (ye + ys)/2.0
            length = mag([ xe - xs, ye - ys ])
            angle = atan(ye - ys, xe - xs)
            loc = angle <= π ? "upper" : "lower"

            influence(x, y, dxdz, dydz) = 1. /(2π*((x - x_14)^2 + (y - y_14)^2))*((y - y_14)*dxdz - (x - x_14)*dydz)
            potential(x, y) = strength/(2π)*atan(y - y_14, x - x_14)
            velocity(x, y) = strength/(2π*((x - x_14)^2 + (y - y_14)^2)).*(y - y_14, -(x - x_14))
            new(xs, ys, xe, ye, x_14, y_14, x_34, y_34, length, angle, loc, strength, vt, cp, influence, potential, velocity)
        end 
    end

    struct VortexPanelSolver2D <: Solution
        panels :: Array{VortexPanel2D, 1}
        uniform :: Uniform2D
        error
        velocity
        potential
        solveStrengths
        tangentialVelocities
        liftCoefficient
        function VortexPanelSolver2D(panels, uniform)
            function solveStrengths()
                """
                Solves for the source strengths of all panels.
                """
                # Construct matrix for normal direction.
                An = [ panel_i.influence(panel_j.x_34, panel_j.y_34, sin(panel_j.angle), cos(panel_j.angle)) for (j, panel_j) in enumerate(panels), (i, panel_i) in enumerate(panels) ]
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
            
            function tangentialVelocities(pressure = true)
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
            new(panels, uniform, error, velocity, potential, solveStrengths, tangentialVelocities, liftCoefficient)
        end
    end

    function panelSolver2D(aeroproblem :: VortexPanelSolver2D, uniform :: Uniform2D)
        strengths = aeroproblem.solveStrengths()
        vts = aeroproblem.tangentialVelocities()
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
        strength
        vt
        cp
        influence
        potential
        velocity
        function DoubletPanel2D(xs, ys, xe, ye, strength=0.0, vt=0.0, cp=0.0)
            xc, yc = (xe + xs)/2.0, (ye + ys)/2.0
            length = mag([ xe - xs, ye - ys ])
            angle = atan(ye - ys, xe - xs)
            loc = angle <= π ? "upper" : "lower"

            influence(x, y, dxdz, dydz) = 1/(2π).*((y/((x - xs)^2 + y^2) - y/((x - xe)^2 + y^2))*dxdz + (-(x - xs)/((x - xs)^2 + y^2) + (x - xe)/((x - xe)^2 + y^2))*dydz)
            potential(x, y) = strength/(4π)*(x - xc)/((x - xc)^2 + (y - yc)^2)
            velocity(x, y) = strength/(2π).*(y/((x - xs)^2 + y^2) - y/((x - xe)^2 + y^2), -(x - xs)/((x - xs)^2 + y^2) + (x - xe)/((x - xe)^2 + y^2))            
            new(xs, ys, xe, ye, xc, yc, length, angle, loc, strength, vt, cp, influence, potential, velocity)
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
        tangentialVelocities
        liftCoefficient
        function DoubletPanelSolver2D(panels, uniform)
            num_panels = length(panels)
            woke_panel = DoubletPanel2D(panels[end].xe, panels[end].ye, 10*panels[end].xe, panels[end].ye)
            function solveStrengths()
                """
                Solves for the source strengths of all panels.
                """
                # Construct doublet matrix for normal direction.
                doublet = [ i == j ? -2/(π*panel_i.length) : panel_i.influence(panel_j.xc, panel_j.yc, -sin(panel_i.angle - panel_j.angle), cos(panel_i.angle - panel_j.angle)) for (j, panel_j) in enumerate(panels), (i, panel_i) in enumerate(panels) ]

                # Construct woke vector
                woke = [ woke_panel.influence(panel_j.xc, panel_j.yc, sin(woke_panel.angle - panel_j.angle), cos(woke_panel.angle - panel_j.angle)) for panel_j in panels ]

                # Kutta condition
                kutta = zeros(num_panels + 1)
                kutta[1] = 1.
                kutta[end-1] = -1.
                kutta[end] = 1.

                An = zeros(num_panels + 1, num_panels + 1)
                An[1:end-1, 1:end-1] = doublet
                An[1:end-1, end] = woke
                An[end, :] = kutta
                
                pans = push!(panels, woke_panel)
                # Construct freestream RHS.
                bn = @. -uniform.magnitude*sin([ uniform.angle + panel.angle for panel in pans ]) 
                
                # Solve system.
                strengths = An\bn
        
                # Update panel strengths.
                for (panel, strength) in zip(pans, strengths)
                    panel.strength = strength
                end

                return strengths
            end
            
            function tangentialVelocities(pressure = true)
                """
                Solves for the velocities of all panels and their pressure coefficients.
                """
                pans = push!(panels, woke_panel)
                # Construct matrix for tangential direction.
                At = [ i == j ? 0. : panel_i.influence(panel_j.xc, panel_j.yc, cos(panel_j.angle), -sin(panel_j.angle)) for (j, panel_j) in enumerate(pans), (i, panel_i) in enumerate(pans) ]

                # Construct freestream RHS.
                bt = @. uniform.magnitude*cos([ uniform.angle + panel.angle for panel in pans ]) 

                # Solve system.
                vts = At*[ panel.strength for panel in pans ] + bt
        
                # Update panel velocities and pressure coefficients.
                for (panel, vt) in zip(pans, vts)
                    panel.vt = vt
                    panel.cp = pressure ? pressureCoefficient2D(0., vt, uniform.magnitude) : 0
                end
        
                return vts
            end
            function liftCoefficient()
                pans = push!(panels, woke_panel)
                strengths = [ panel.strength for panel in pans ]
                cl = (1/uniform.magnitude*sum(strengths[1:end-2] - strengths[2:end-1]), -strengths[end]/uniform.magnitude)
                return cl
            end
            pans = push!(panels, woke_panel)
            # Compute errors
            error() = sum([ panel.strength*panel.length for panel in panels ])
            # Compute velocities
            velocity(x, y) = uniform.velocity(x,y) .+ foldl( (u, v) -> u .+ v, [ panel.velocity(panel.strength, x, y) for panel in pans ] )
            # Compute potential
            potential(x, y) = uniform.potential(x, y) .+ sum([ panel.potential(x, y) for panel in panels ])
            new(panels, uniform, error, velocity, potential, solveStrengths, tangentialVelocities, liftCoefficient)
        end
    end

    function panelSolver2D(aeroproblem :: DoubletPanelSolver2D, uniform :: Uniform2D)
        strengths = aeroproblem.solveStrengths()
        vts = aeroproblem.tangentialVelocities()
        cl = aeroproblem.liftCoefficient()
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
        tangentialVelocities
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
                bn = push!(@. -uniform.magnitude * cos([ uniform.angle - panel.angle for panel in panels ]), -uniform.magnitude * (sin(uniform.angle - panels[1].angle) + sin(uniform.angle - panels[end].angle)))

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
            function tangentialVelocities(pressure=true)
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
            new(panels, uniform, gamma, error, velocity, potential, solveStrengths, tangentialVelocities, liftCoefficient)
        end
    end

mutable struct DoubletSourcePanel2D <: Solution
    xs
    ys
    xe 
    ye
    xc
    yc
    xsl
    ysl
    xel
    yel
    xcl
    ycl
    length
    angle
    loc :: T where T <: AbstractString
    strength
    vt
    cp
    sourceInfluence
    doubletInfluence
    potential
    velocity
    function DoubletSourcePanel2D(xs, ys, xe, ye, strength=0.0, vt=0.0, cp=0.0)
        # Properties
        length = mag([xe - xs, ye - ys])
        angle = atan(ye - ys, xe - xs)
        loc = angle <= π ? "upper" : "lower"
        xc, yc = (xe + xs)/2.0, (ye + ys)/2.0
        # Global coordinates
        # if 0 < angle <= π/2
        #     xc, yc = (xe + xs)/2.0 + 0.05, (ye + ys)/2.0 - 0.05
        # elseif π/2 < angle <= π
        #     xc, yc = (xe + xs)/2.0 - 0.05, (ye + ys)/2.0 - 0.05
        # elseif -π < angle <= -π/2
        #     xc, yc = (xe + xs)/2.0 - 0.05, (ye + ys)/2.0 + 0.05
        # # elseif -π/2 < angle <= 0
        # else
        #     xc, yc = (xe + xs)/2.0 + 0.05, (ye + ys)/2.0 + 0.05
        # end

        # Panel coordinates
        xsl, ysl, xel, yel = 0., 0., length, 0.
        xcl, ycl = length/2., 0

        # Influence coefficients, defined in panel coordinates
        sourceInfluence(x, y) = 1/(4π)*((x - xsl)log((x - xsl)^2 + y^2)) - (x - xel)log((x - xel)^2 + y^2) + y*(atan(y, x - xel) - atan(y, (x - xsl)))
        doubletInfluence(x, y) = -1/(2π)*(atan(y, x - xel) - atan(y, x - xsl))

        # Fluid methods (IGNORE)
        potential(x, y) = strength/(4π)*(x - xc)/((x - xc)^2 + (y - yc)^2)
        velocity(x, y) = strength/(2π).*(y/((x - xs)^2 + y^2) - y/((x - xe)^2 + y^2), -(x - xs)/((x - xs)^2 + y^2) + (x - xe)/((x - xe)^2 + y^2))            
        new(xs, ys, xe, ye, xc, yc, xsl, ysl, xel, yel, xcl, ycl, length, angle, loc, strength, vt, cp, sourceInfluence, doubletInfluence, potential, velocity)
    end 
end

struct DoubletSourcePanelSolver2D <: Solution
    """
    Applies the vortex source panel method to a list of panels and a uniform flow.
    """
    panels :: Array{DoubletSourcePanel2D, 1}
    uniform :: Uniform2D
    gamma
    error
    velocity
    potential
    solveStrengths
    tangentialVelocities
    liftCoefficient
    function DoubletSourcePanelSolver2D(panels, uniform, gamma=0.0)
        num_panels = length(panels)
        woke_panel = DoubletSourcePanel2D(panels[end].xe, panels[end].ye, 1000*panels[end].xe, panels[end].ye)
        print(woke_panel.xs, woke_panel.ys, woke_panel.xe, woke_panel.ye)

        function solveStrengths()
            """
            Solves for the strengths of all panels.
            """
            # Wake vector
            woke_vector = [ woke_panel.doubletInfluence(panelCoords(panel_i.xc, panel_i.yc, woke_panel.xs, woke_panel.ys, woke_panel.angle)...) for panel_i in panels ]
            
            # Doublet matrix
            doublet_matrix = [ i == j ? 0.5 : panel_j.doubletInfluence(panelCoords(panel_i.xc, panel_i.yc, panel_j.xs, panel_j.ys, panel_j.angle)...) for (j, panel_j) in enumerate(panels), (i, panel_i) in enumerate(panels) ]
            print(doublet_matrix, '\n')
            
            # LHS
            An = copy(doublet_matrix)
            An[:, 1] = doublet_matrix[:, 1] - woke_vector
            An[:, end] = doublet_matrix[:, end] + woke_vector 
            print(An, '\n')

            # Source matrix
            source_matrix = [ i == j ? 0.5*panel_j.sourceInfluence(panelCoords(panel_i.xc, panel_i.yc, panel_j.xs, panel_j.ys, panel_j.angle)...) : panel_j.sourceInfluence(panelCoords(panel_i.xc, panel_i.yc, panel_j.xs, panel_j.ys, panel_j.angle)...) for (j, panel_j) in enumerate(panels), (i, panel_i) in enumerate(panels) ]
            # print(source_matrix, '\n')
        
            # Source vector
            source_vector = -@. uniform.magnitude*sin([ uniform.angle - panel.angle for panel in panels ]) 
            # print(source_vector, '\n')
            
            # RHS
            b = -source_matrix*source_vector
            # print(b)

            # Solve system
            strengths = An\b

            # Update panel strengths.
            for (panel, strength) in zip(panels, strengths)
                panel.strength = strength
            end

            return strengths
        end
        function tangentialVelocities(pressure=true)
            """
            Solves for the velocities of all panels and their pressure coefficients.
            """
            for panel in panels
                panel.strength = panel.xc*cos(uniform.angle) + panel.yc*sin(uniform.angle) + panel.strength
            end
            vts = [ (panel1.strength - panel2.strength)/mag([panel1.xc - panel2.xc, panel1.yc - panel2.yc]) + uniform.magnitude for (panel1, panel2) in zip(panels[1:end-1], panels[2:end])]
            # Update panel velocities and pressure coefficients.
            for (panel, vt) in zip(panels, vts)
                panel.vt = vt
                panel.cp = pressure ? pressureCoefficient2D(0., vt, uniform.magnitude) : 0.
            end

            return vts
        end
        function liftCoefficient()
            # Lift coefficient
            xs = [ panel.xs for panel in panels ]
            c = abs(maximum(xs) - minimum(xs))
            cl = sum([ -panel.cp*panel.length*cos(panel.angle) for panel in panels ])/c

            return cl
        end
        # Compute errors
        error() = sum([ panel.strength*panel.length for panel in panels ])
        # Compute velocities
        velocity(x, y) = uniform.velocity(x,y) .+ foldl( (u, v) -> u .+ v, [ panel.velocity(panel.strength, x, y) for panel in panels ] )
        # Compute potential
        potential(x, y) = uniform.potential(x, y) .+ sum([ panel.potential(x, y) for panel in panels ])
        new(panels, uniform, gamma, error, velocity, potential, solveStrengths, tangentialVelocities, liftCoefficient)
    end
end

    function panelSolver2D(aeroproblem :: DoubletSourcePanelSolver2D, uniform :: Uniform2D)
        strengths = aeroproblem.solveStrengths()
        vts = aeroproblem.tangentialVelocities()
        cl = aeroproblem.liftCoefficient()
        error = aeroproblem.error()
        return strengths, vts, cl, error
    end

    function gridData(objects :: Array{<:Solution, 1}, xs)
        vels = foldl( (v1,v2) -> [ u .+ v for (u,v) in zip(v1, v2)], [ velocity(object, xs) for object in objects ])
        pots = foldl( (v1,v2) -> v1 .+ v2, [ potential(object, xs) for object in objects ])

        return vels, pots
    end

    function gridData(object :: Solution, xs)
        vels = velocity(object, xs)
        pots = potential(object, xs)
        
        return vels, pots
    end 
    
    # Transforms (x, y) to the coordinate system with (x_s, y_s) as origin oriented at angle_s.
    panelCoords(x, y, x_s, y_s, angle_s) = rotation(x - x_s, y -  y_s, angle_s)

    # Rotation matrices
    invRotation(x, y, angle) = (x*cos(angle) - y*sin(angle), x*sin(angle) + y*cos(angle))
    rotation(x, y, angle) = (x*cos(angle) + y*sin(angle), -x*sin(angle) + y*cos(angle))

    pressureCoefficient2D(vels :: Array, freestream) = pressureCoefficient2D.(map(first, vels), map(last, vels), freestream)

    velocity(object :: Solution, xs) = [ object.velocity(x...) for x in xs ] 
    potential(object :: Solution, xs) = [ object.potential(x...) for x in xs ]
    
    # Euclidean norm of a vector
    mag(x) = sqrt(sum(s -> s^2, x))

    # Infinite vortices
    infiniteVortices(strength :: Float64, spacing :: Float64, x :: Float64, y :: Float64) = (strength/(2*spacing)*sinh(2π*y/spacing)/(cosh(2π*y/spacing) - cos(2π*x/spacing)), strength/(2*spacing)*sin(2π*x/spacing)/(cosh(2π*y/spacing) - cos(2π*x/spacing)))

    # Joukowski transformation
    joukowski(z :: Complex, c) = z + c^2/z

    # Pressure coefficient
    pressureCoefficient2D(u, v, freestream_speed) = 1. - (u^2 + v^2)/freestream_speed^2

    function cosinePanels(x :: Array{<:Real}, y :: Array{<:Real}, n :: Integer = 40)
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

        panels = [ SourcePanel2D(x_ends[i], y_ends[i], x_ends[i+1], y_ends[i+1]) for i in 1:n ]

        return panels
    end
end