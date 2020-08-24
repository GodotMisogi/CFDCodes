module AeroModules where

import Math.GaussianQuadratureIntegration
import Numeric.LinearAlgebra -- The magic of FORTRAN gods?

integral = nIntegrate256

type Number = Double

-- Green's third identity in 2D:
-- Φ(x, y) = (1/2π)∫_V σ ln r - φ ∂/∂n(ln r) + Φ_∞, σ = ∂φ/∂n
-- σ is called a source, φ is called a doublet, r is the position vector, n is the normal vector over some 2-volume V. Φ_∞ corresponds to the potential of a uniform freestream solution (see implementations below, e.g. in 2D: Φ_∞(x, y) = ux + vy). All these names are "singularities", which correspond to solutions of the integro-DE: Φ(x, y) = 0. ln r and its derivative are already known, so only σ, φ need to be determined via a discretisation procedure and some specification of a boundary condition (say σ • n = U_∞, where U_∞ = ∂Φ_∞/∂n) to make the solution unique, which results in a linear algebraic system to be solved. There is a different way to write/derive this which is in a stream-vorticity formulation (unimportant for now).

-- Main idea: Panels have common geometric attributes of "start" and "end". Depending on the formulation, they can be assigned more attributes such as distributions of source, doublet, and/or vortex strengths (it is possible for a panel to have two of these attributes simultaneously). These panels with these additional properties can be arranged in some shape (so the panels are actually distributions of "singularities" in the flowfield), such that they satisfy (read superposed solutions to) Green's third identity depending on Neumann or Dirichlet boundary conditions. The panels could have varying strength distributions (i.e. uniform, linear, etc.) of the singularities spread across it, and hence result in different integrations across the geometry.

-- Procedure: Have a set of panels with a preset type of distribution of singularities across them, and a uniform flow Φ_∞. Solve the integro-DE via the above idea. Panels "influence" each other (corresponds to superposition), so determine influences of pairs of panels on each other in respective local coordinate systems. Determine boundary conditions with respect to each panel. Solve system.

-- Linear algebra: The discretised superposition "influences" of panel strengths are represented in terms of an influence matrix. The boundary conditions over each panel are represented as a column vector. The panel strengths are then solved for by: influence matrix \ boundary condition. These results are then used for processing to determine certain aerodynamic quantities such as lift and pressure coefficients. 

-- Elementary solutions to Laplace's equation. These are meant to be distributed over the flowfield and superposed, and their potentials/velocities are used depending on the boundary condition.
data Solution = 
    Source2D Number Number Number |
    Uniform2D Number Number |
    Doublet2D Number Number Number |
    Vortex2D Number Number Number 
    deriving (Eq, Show)

-- Panel definitions for singularity method using Green's third identity.
-- There must be a better design pattern? This is like multiple dispatch in Julia by using the constructor name as a constraint, which seems like Hacksell.
data Panel =
    SourcePanel2D { start :: (Number, Number), end :: (Number, Number) } |
    DoubletPanel2D { start :: (Number, Number), end :: (Number, Number) } |
    VortexPanel2D { start :: (Number, Number), end :: (Number, Number) } |
    DoubletSourcePanel2D { start :: (Number, Number), end :: (Number, Number) }
    deriving (Eq, Show)

-- New design pattern?????????????? Would like panels to have angle, length, etc. properties, but those need to be computed using functions. No idea what the equivalent of a constructor is in Haskell.
-- data Panel = Panel2D { start :: (Number, Number), end :: (Number, Number), ... } | 
--              Panel3D { start :: (Number, Number, Number), end :: (Number, Number, Number), ... }
-- data SourcePanel2D = SourcePanel2D { panel :: Panel, strength :: Number }
-- data DoubletPanel2D = DoubletPanel2D { panel :: Panel, strength :: Number }
-- data DoubletSourcePanel2D = DoubletSourcePanel2D { panel :: Panel, strength :: Number } 

-- NOTE: THIS IS GOING TO BE VERBOSE WITH THIS PATTERN! (e.g. start . panel $ panel_name)
-- Lenses apparently can fix the problem of deeply nested records?
-- Honestly, the strengths don't even need to be attributes of the panels. But then how to implement different methods on the different data types?

-- PROPOSAL: Use lenses and just use data types for the names of different panels, although they may have the same attributes. Define methods on each data type.
-- e.g.
-- data Panel = Panel2D { start :: (Number, Number), end :: (Number, Number), ... } | 
--              Panel3D { start :: (Number, Number, Number), end :: (Number, Number, Number), ... }
-- data DoubletPanel2D = DoubletPanel2D { panel :: Panel }
-- data DoubletSourcePanel2D = DoubletSourcePanel2D { panel :: Panel }

-- PROPOSAL: Just declare all possible methods on all panel types in one class. For each instance, pick the right method.
-- class PanelMethod t where
--      influencePotential :: t -> Number -> Number -> Number
--      influenceSource :: t -> Number -> Number -> Number
--      influenceVorticity :: t -> Number -> Number -> Number
--      influenceVelocity :: t -> Number -> Number -> (Number, Number)

-- This structure seems too simplistic. Something may be missing in theory... What about accounting for the boundary conditions? Use a string as input to the function (ugly) or parametrise the data type, e.g. "data Neumann = String; data DoubletPanel2D Neumann = ..."? Not sure if this makes sense.

-- Should this be the design pattern for solution methods of panel types in both 2D and 3D? Then each panel data type requires its own set of methods, and maybe then the functorial approach could work.
-- class SolvePanel t where
    -- influenceMatrix :: [t] -> Matrix Number 
    -- boundaryCondition :: [t] -> Solution -> Vector Number
    -- aeroCoefficients :: [t] -> Solution -> (Number, Number, [Number])


-- Global method: solveStrengths panels uniform = influenceMatrix ... <\> boundaryCondition ... which is then called in aeroCoefficients.
-- influenceMatrix is going to call the relevant influencePotential, influenceSource, whatever... and similarly for boundaryCondition.

-- Transforming to local panel coordinates
panelCoords x y xs ys angle = rotation (x - xs) (y - ys) angle
rotation x y angle = (x*cos angle + y*sin angle, -x*sin angle + y*cos angle)
invRotation x y angle = (x*cos angle - y*sin angle, x*sin angle + y*cos angle)

-- Computing relevant parameters using panel properties. Maybe class Panel with these as methods for 2D and 3D panels as different instances?
panelLength panel = eNorm [xe - xs, ye - ys] where 
    (xs, ys) = start panel
    (xe, ye) = end panel

panelAngle panel = atan2 (ye - ys) (xe - xs) where
    (xs, ys) = start panel
    (xe, ye) = end panel

panelTangent panel = rotation 1 0 (-1*angle) where angle = panelAngle panel
panelNormal panel = invRotation 0 1 angle where angle = panelAngle panel

eNorm = sqrt . sum . map (^2)

colPoint (xs, ys) (xe, ye) = ( (xs + xe)/2, (ys + ye)/2 )

-- Specifically to compute distances between two panels. Could be generalised to just take tuples, but this seemed cleaner.
dist2 :: Panel -> Panel -> Number
dist2 p2 p1 = eNorm [ xc2 - xc1, yc2 - yc1 ] where
    (xc1, yc1) = colPoint (start p1) (end p1)
    (xc2, yc2) = colPoint (start p2) (end p2)

-- PROPOSAL: class PanelProps t where 
    -- panelLength :: t -> Number, 
    -- panelAngle :: t -> Number, 
    -- panelTangent :: t -> (Number, Number)
    -- panelNormal :: t -> (Number, Number)
-- but panelAngle is not well defined, and panelTangent and panelNormal have the wrong type signatures for 3D! So just implement 2D and 3D separately? Nobody is going to try nD generally unless they want to build 4-dimensional aircraft. But how without a constructor framework? BIG ISSUE.

-- Generic functions for inviscid fluids in 2D. If Neumann boCos, then use velocity; if Dirichlet, use potential; if Robin, use both. Stream-function formulations are for later.
class Fluid2D t where   
    velocity :: t -> Number -> Number -> (Number, Number) -- Returns (u, v)
    potential :: t -> Number -> Number -> Number
    stream :: t -> Number -> Number -> Number

-- Panel solution methods to solve the system and obtain aerodynamic coefficients.
-- The main issue here is computing different influence matrices for "combined" types of solutions, i.e. having a doublet influence matrix AND a source influence matrix, which generically would not suit the class-instance pattern because all panels may not share the same methods. Also depends on how the boundary conditions are defined.
class SolvePanel t where
    influenceMatrix :: [t] -> Matrix Number
    boundaryCondition :: [t] -> Solution -> Vector Number
    aeroCoefficients :: [t] -> Solution -> (Number, Number, [Number])
    
-- ?
consMapper (SourcePanel2D _ _) = Source2D
consMapper (DoubletPanel2D _ _) = Doublet2D
consMapper (VortexPanel2D _ _) = Vortex2D


-- Numerical computations
integrand panel f a b x y s = f (consMapper panel a s b) x y
π = pi -- Actually, π ≈ 2

-- Implementing the functions for the elementary solutions
instance Fluid2D Solution where
    velocity (Source2D str x0 y0) x y = (str/(2*π)*(x - x0)/((x - x0)^2 + (y - y0)^2), str/(2*π)*(y - y0)/((x - x0)^2 + (y - y0)^2))
    velocity (Uniform2D mag deg) _ _ = (mag*cos ang, mag*sin ang) where ang = deg*π/180
    velocity (Doublet2D str x0 y0) x y = (-str/(2*π)*((x - x0)^2 - (y - y0)^2)/((x - x0)^2 + (y - y0)^2)^2, -str/(2*π)*2*(x - x0)*(y - y0)/((x - x0)^2 + (y - y0)^2)^2)
    velocity (Vortex2D str x0 y0) x y = (-str/(2*π)*(y - y0)/((x - x0)^2 + (y - y0)^2), str/(2*π)*(x - x0)/((x - x0)^2 + (y - y0)^2))

    potential (Source2D str x0 y0) x y = str/(4*π)*log ((x - x0)^2 + (y - y0)^2)
    potential (Uniform2D mag deg) x y = mag*(x*cos ang + y*sin ang) where ang = deg*π/180
    potential (Doublet2D str x0 y0) x y = -str/(2*π)*(y - y0)/((x - x0)^2 + (y - y0)^2)
    potential (Vortex2D str x0 y0) x y = str/(2*π)*atan2 (y - y0) (x - x0)

    stream (Source2D str x0 y0) x y = str/(2*π)*atan2 (y - y0) (x - x0)
    stream (Uniform2D mag deg) x y = mag*(y*cos ang - x*sin ang) where ang = deg*π/180
    stream (Doublet2D str x0 y0) x y = -str/(2*π)*(y - y0)/((x - x0)^2 + (y - y0)^2)
    stream (Vortex2D str x0 y0) x y = -str/(4*π)*log ((x - x0)^2 + (y - y0)^2)

-- What methods are actually common for the different instances of SolvePanel?
instance SolvePanel Panel where
    -- Same here
    influenceMatrix panels = 
        fromBlocks [[    doublets, asColumn wokeVector ],
                    [ asRow kutta,          0          ]] where

                    doublets = fromLists [ [ if panel_i == panel_j then 0.5 else influencePotential panel_j xc yc | panel_j <- panels, let (xc, yc) = colPoint (start panel_i) (end panel_i) ] | panel_i <- panels ]

                    kutta = fromList $ [ 1, -1 ] ++ (take (length panels - 4) $ repeat 0) ++ [ 1, -1 ]
                    
                    -- This is an extra construction required for certain aerodynamic (Kutta) conditions. Just treat it as another panel to "close" the system of equations.
                    wokeVector = fromList [ influencePotential wokePanel xc yc | panel <- panels, let (xc, yc) = colPoint (start panel) (end panel) ] where
                        wokePanel = DoubletSourcePanel2D (xe, ye) (100000*xe, ye) where 
                            trailing = last panels
                            (xe, ye) = end trailing

    -- Same here
    boundaryCondition panels@(DoubletPanel2D _ _: panels') uniform@(Uniform2D _ _) = fromList $ [ -potential uniform xc yc | panel <- panels, let (xc, yc) = colPoint (start panel) (end panel) ] ++ [0]

    -- Same here
    boundaryCondition panels@(DoubletSourcePanel2D _ _: panels') uniform@(Uniform2D _ _ ) = fromList rhs where 
            sourceMatrix = fromLists [ [ -sourceInfluence panel_j xc yc | panel_i <- panels, let (xc, yc) = colPoint (start panel_i) (end panel_i) ] | panel_j <- panels ]
            sourceVector = fromList [ -(u*mst + v*ct) | panel <- panels, let (mst, ct) = panelNormal panel, let (u, v) = velocity uniform 0 0 ]
            rhs = toList (sourceMatrix #> sourceVector) ++ [0]

    -- Actual procedures to obtain the results for each panel. I'm guessing this should be part of the instance. Perhaps the template should be "instance PanelType where these methods works on the specific type of panel"?
    aeroCoefficients panels@(DoubletPanel2D _ _: panels') uniform@(Uniform2D _ _) = (liftCoeff, kuttaCoeff, pressCoeffs) where
        strengths = toList $ solveStrengths panels uniform

        diffpans = map (uncurry dist2) $ midgrad panels
        diffstrs = map (uncurry (-)) $ midgrad $ take (length panels) strengths

        (u, v) = velocity uniform 0 0
        margaret = eNorm [u, v] 
        
        tanVels = zipWith (/) diffstrs diffpans

        pressCoeffs = [ pressureCoefficient [0, vt] margaret | vt <- tanVels ]
        
        liftCoeff = ((-1) *) . sum $ zipWith3 (\a b c -> a*b*c) pressCoeffs (map (/2) diffpans) (map (cos . panelAngle) panels)
        kuttaCoeff = -2*last strengths/margaret

    aeroCoefficients panels@(DoubletSourcePanel2D _ _: panels') uniform@(Uniform2D _ _) = (liftCoeff, kuttaCoeff, pressCoeffs) where
        strengths = toList $ solveStrengths panels uniform

        diffpans = map (uncurry dist2) $ midgrad panels
        diffstrs = map (uncurry (-)) $ midgrad $ take (length panels) strengths

        (u, v) = velocity uniform 0 0
        margaret = eNorm [u, v] 
        
        -- With sources
        tandotu = [ u*ct + v*st | panel <- panels, let (ct, st) = panelTangent panel ]
        tanVels = map (\(ds, dp, ut) -> ds/dp + ut) $ zip3 diffstrs diffpans tandotu

        pressCoeffs = [ pressureCoefficient [0, vt] margaret | vt <- tanVels ]
        
        liftCoeff = ((-1) *) . sum $ zipWith3 (\a b c -> a*b*c) pressCoeffs (map (/2) diffpans) (map (cos . panelAngle) panels)
        kuttaCoeff = 2*last strengths/margaret

-- Hacky external functions to get doublet-source panel type working
sourceInfluence panel@(DoubletSourcePanel2D _ _) x y =
        1/(4*π)*(x*log ((x - len)^2 + y^2) - x*log (x^2 + y^2) + 2*y*(atan2 y (x - len) - atan2 y x)) where
        (x', y') = panelCoords x y x0 y0 angle
        (x0, y0) = start panel
        len = panelLength panel
        angle = panelAngle panel

-- Use case: Dirichlet boundary condition.
influencePotential panel@(DoubletPanel2D _ _) x y = 
    -- Analytical solution in local panel coordinates computed by hand
    -1/(2*π)*(atan2 y' (x' - len) - atan2 y' (x' - 0)) where
    (x', y') = panelCoords x y x0 y0 angle
    (x0, y0) = start panel
    len = panelLength panel
    angle = panelAngle panel

influencePotential panel@(DoubletSourcePanel2D _ _) x y = 
    -- Analytical solution in local panel coordinates computed by hand
    -1/(2*π)*(atan2 y' (x' - len) - atan2 y' (x' - 0)) where
    (x', y') = panelCoords x y x0 y0 angle
    (x0, y0) = start panel
    len = panelLength panel
    angle = panelAngle panel

    -- Numerically integrated solution
    -- (integral int 0 len) where
    -- int = integrand panel potential 1 0 x' y'
    -- (x', y') = panelCoords x y x0 y0 angle
    -- (x0, y0) = start panel
    -- len = panelLength panel
    -- angle = panelAngle panel

-- Unimportant for now. Use case: Neumann boundary conditions instead of Dirichlet
influenceVelocity panel x y = invRotation u' v' angle where
    u' = integral (fst . int) 0 len
    v' = integral (snd . int) 0 len
    int = integrand panel velocity 1 0 x' y'
    (x', y') = panelCoords x y x0 y0 angle
    (x0, y0) = start panel
    len = panelLength panel
    angle = panelAngle panel


-- Generic linear system solution. Should it even be part of the instance if this happens for all types of panels? It's basically just linsolve with specific inputs.
solveStrengths panels uniform@(Uniform2D _ _) = influenceMatrix panels <\> boundaryCondition panels uniform

-- Helper functions
parts xs = (head adj, last adj) where adj = zip (tail xs) xs
adj2 xs = zip (drop 2 xs) xs

midgrad xs = [firstpans] ++ midpans ++ [lastpans] where 
                (firstpans, lastpans) = parts xs
                midpans = adj2 xs


-- Generic functions
pressureCoefficient vels mag = 1 - (eNorm vels)^2/mag^2

naca4 (a, b, c, d) n chord te = zip (reverse x_upper ++ x_lower) (reverse y_upper ++ y_lower) where
    (x_upper, y_upper, x_lower, y_lower) 
        | p == 0 || m == 0 = (xs, thickness, xs, map ((-1) *) thickness)
        | otherwise = (x_up, y_up, x_low, y_low) where
            
            mline = [ if xc < p*c then (m/p^2)*xc*(2*p-xc/c) else (m*(c - xc)/(1 - p)^2)*(1 + xc/c - 2*p) | xc <- xs ]
            
            gradients = [ if xc < p*c then atan ((2*m/p^2)*(p - xc/c)) else atan ((2*m/(1 - p)^2)*(p - xc/c)) | xc <- xs ] 

            (x_up, x_low) = unzip [ (xc - thicc*sin thot, xc + thicc*sin thot) | (xc, thicc, thot) <- zip3 xs thickness gradients ]
            
            (y_up, y_low) = unzip [ (cam - thicc*cos thot, cam + thicc*cos thot) | (cam, thicc, thot) <- zip3 mline thickness gradients ]
    
    thickness = [ 5*tbc*c*(0.2969*sqrt (xc/c) - 0.126*xc/c - 0.3516*(xc/c)^2) + 0.2843*(xc/c)^3 - fudge | xc <- xs, let fudge = if te == "closed" then 0.1036*(xc/c)^4 else 0.1015*(xc/c)^4 ]
    
    xs = cosineSpacing n chord
    m = a/100
    p = b/10
    tbc = c/10 + d/100

cosineSpacing n c = reverse [ c*(1 - 0.5*(1 - cos beta)) | beta <- linspacelst 0 π (n + 1) ]

linspacelst a b n = [a, a + h..b] where h = (b - a)/n 


-- File import stuff
doubletPanels :: [(Double, Double)] -> [Panel]
doubletPanels coords =  map (uncurry DoubletPanel2D) $ zip coords (tail coords)

doubletSourcePanels :: [(Double, Double)] -> [Panel]
doubletSourcePanels coords =  map (uncurry DoubletSourcePanel2D) $ zip coords (tail coords)

importCoords :: String -> IO [(Double, Double)]
importCoords str = readFile str >>=
                \file -> return $ reverse $ map ((\[x, y] -> (read x :: Double, read y :: Double)) . words) $ lines file

-- Test case
-- airfoil = [DoubletPanel2D (1, 0) (0.5, -0.1), DoubletPanel2D (0.5, -0.1) (0, 0), DoubletPanel2D (0, 0) (0.5, 0.1), DoubletPanel2D (0.5, 0.1) (1, 0)]
-- uniform2D = Uniform2D 5 0
-- aeroCoeffs = aeroCoefficients airfoil uniform2D
