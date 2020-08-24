module panelSolver

abstract type Solution end

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
    panels :: Array{DoubletSourcePanel2D, 1}
    uniform :: Uniform2D
    error
    velocity
    potential
    solveStrengths
    aerodynamicsss
    liftCoefficient
    function DoubletSourcePanelSolver2D(panels, uniform, sources=true, vel_kutta=true)
        """
        Solves Laplace's equation over an airfoil using constant-strength source and doublets based on a Dirichlet boundary condition with the Morino-Kutta condition.

        Note: The method can be reduced to the constant-strength doublet potential method by setting the third argument to 'false'.
        """
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
            
            kutta = zeros(num_panels + 1)
            if vel_kutta
                # Morino's velocity Kutta condition
                kutta[1] = 1.
                kutta[2] = -1.
                kutta[end-2] = 1.
                kutta[end-1] = -1.
            else
                # Woke Kutta condition
                kutta[1] = 1.
                kutta[end-1] = -1.
                kutta[end] = 1.
            end
            # LHS with explicit Kutta condition
            An = zeros(num_panels + 1, num_panels + 1)
            An[1:end-1, 1:end-1] = doublet_matrix
            An[1:end-1, end] = woke_vector 
            An[end, :] = kutta
            println(An)

            # LHS with implicit Kutta condition
            # An = zeros(num_panels, num_panels)
            # An = doublet_matrix
            # An[:, 1] = doublet_matrix[:, 1] - woke_vector 
            # An[:, end] = doublet_matrix[:, end] + woke_vector

            b = zeros(num_panels + 1)
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
                b[1:end-1] = -source_matrix*source_vector
            else
                # Freestream RHS
                b[1:end-1] = [ -sum(u.*(panel.xc, panel.yc)) for panel in panels ]
            end
            println(b)

            # Solve system
            strengths = An\b
            # println(strengths)

            # Update panel doublet strengths
            for (panel, strength) in zip(panels, strengths)
                panel.doublet_strength = strength
            end

            woke_panel.doublet_strength = strengths[end]

            return strengths
        end
        function aerodynamicsss()
            """
            Solves for the velocities of all panels and their pressure coefficients.
            """
            
            xs = [ panel.xs for panel in panels ]
            c = abs(maximum(xs) - minimum(xs))

            # Update panel velocities, pressure and lift coefficients.
            cl = 0
            for (i, panel) in enumerate(panels)
                if panel == panels[1]
                    R = mag([panels[2].xc - panels[1].xc, panels[2].yc - panels[1].yc])
                    panel.vt = sources ? (panels[2].doublet_strength - panels[1].doublet_strength)/R + sum(u.*panel.tangent) : (panels[2].doublet_strength - panels[1].doublet_strength)/R
                elseif panel == panels[end]
                    R = mag([panels[end].xc - panels[end-1].xc, panels[end].yc - panels[end-1].yc])
                    panel.vt = sources ? (panels[end].doublet_strength - panels[end-1].doublet_strength)/R + sum(u.*panel.tangent) : (panels[end].doublet_strength - panels[end-1].doublet_strength)/R 
                else
                    R = mag([panels[i+1].xc - panels[i-1].xc, panels[i+1].yc - panels[i-1].yc])
                    panel.vt = sources ? (panels[i+1].doublet_strength - panels[i-1].doublet_strength)/R + sum(u.*panel.tangent) : (panels[i+1].doublet_strength - panels[i-1].doublet_strength)/R
                end
                panel.cp = pressureCoefficient2D(0., panel.vt, uniform.magnitude)
                cl -= panel.cp*(R/2)*cos(panel.angle)/c
            end
            cl_woke = 2*woke_panel.doublet_strength/uniform.magnitude
            cl =  (cl, sources ? cl_woke : -cl_woke )
            vts = [ panel.vt for panel in panels ]
            cps = [ panel.cp for panel in panels ]
            # println(cps)

            return (cps, cl)
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

function panelSolver2D(panels :: Array{DoubletSourcePanel2D, 1}, uniform :: Uniform2D, sources = true, vel_kutta=true)
    airfoil = DoubletSourcePanelSolver2D(panels, uniform, sources, vel_kutta)
    strengths = airfoil.solveStrengths()
    vts, cl = airfoil.aerodynamicsss()
    error = airfoil.error()
    return airfoil, strengths, vts, cl, error
end

# Performs velocity and potential calculations on a grid
# function gridData(objects :: Array{<:Solution, 1}, xs, pressure = true)
#     vels = foldl( (v1, v2) -> [ u .+ v for (u, v) in zip(v1, v2)], [ velocity(object, xs) for object in objects ])
#     pots = foldl( (v1, v2) -> v1 .+ v2, [ potential(object, xs) for object in objects ])
# end

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

# Cavity termination models 
cavPressure(sf, sl, A = 0.5, ν = 1.0, λ = 0.1) = sf < (1 - λ)*sl ? 0 : A*((sf - (1 - λ)*sl)/(sl - (1 - λ)*sl))^ν

# Euclidean norm of a vector
mag(x) = sqrt(sum(s -> s^2, x))

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
        while j < length(x) - 1
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

    # x, y = push!(x, x[1]), push!(y, y[1])

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

function NACA4(digits :: Tuple, c, n, closed_te=false, split=false)
    
    # Airfoil characteristics
    # Camber
    m = digits[1]/100
    # Position
    p = digits[2]/10
    # Thickness-to-chord ratio
    t_by_c = (10*digits[3] + digits[4])/100

    # Cosine spacing
    angles = range(0, stop = π, length = n + 1)
    xs = reverse([ c*(1 - 0.5*(1 - cos(beta))) for beta in angles ] )

    # Thickness distribution
    thickness = closed_te ? [ 5*t_by_c*c*(0.2969*sqrt(xc/c) - 0.126*xc/c - 0.3516*(xc/c)^2 + 0.2843*(xc/c)^3 - 0.1036*(xc/c)^4) for xc in xs ] : [ 5*t_by_c*c*(0.2969*sqrt(xc/c) - 0.126*xc/c - 0.3516*(xc/c)^2 + 0.2843*(xc/c)^3 - 0.1015*(xc/c)^4) for xc in xs ] 
    
    if p == 0.0 || m == 0.0
        x_upper = xs
        y_upper = thickness
        x_lower = xs
        y_lower = -thickness
    else
        # Compute m
        mline = [ xc < p*c ? (m/p^2)*xc*(2*p-xc/c) : (m*(c - xc)/(1 - p)^2)*(1 + xc/c - 2*p) for xc in xs ]
        # Compute gradients
        gradients = [ xc < p*c ? atan((2*m/p^2)*(p - xc/c)) : atan((2*m/(1 - p)^2)*(p - xc/c)) for xc in xs ] 
        # Upper surface
        x_upper = [ xc - thicc*sin(thot) for (xc, thicc, thot) in zip(xs, thickness, gradients) ] 
        y_upper = [ cam + thicc*cos(thot) for (cam, thicc, thot) in zip(mline, thickness, gradients) ]
        # Lower surface
        x_lower = [ xc + thicc*sin(thot) for (xc, thicc, thot) in zip(xs, thickness, gradients) ] 
        y_lower = [ cam - thicc*cos(thot) for (cam, thicc, thot) in zip(mline, thickness, gradients) ]
    end

    if split
        upper = [ (x, y) for (x, y) in zip(x_upper, y_upper) ]
        lower = [ (x, y) for (x, y) in zip(x_lower, y_lower) ]
        return (reverse(upper), lower)
    else  
        (X, Y) = append!(reverse(x_upper), x_lower), append!(reverse(y_upper), y_lower)
        return (X, Y)
    end
end

end