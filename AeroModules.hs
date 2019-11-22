import Math.GaussianQuadratureIntegration
import Numeric.LinearAlgebra
import Data.Char

integral = nIntegrate256

type Number = Float

data Solution = 
    Source2D Number Number Number |
    Uniform2D Number Number |
    Doublet2D Number Number Number |
    Vortex2D Number Number Number 
    deriving (Eq, Show)

data Panel =
    SourcePanel2D { start :: (Number, Number), end :: (Number, Number) } |
    DoubletPanel2D { start :: (Number, Number), end :: (Number, Number) } |
    VortexPanel2D { start :: (Number, Number), end :: (Number, Number) }
    deriving (Eq, Show)

class Fluid2D t where
    velocity :: t -> Number -> Number -> (Number, Number)
    potential :: t -> Number -> Number -> Number
    stream :: t -> Number -> Number -> Number

class SolvePanel t where
    solveStrengths :: [t] -> Vector Number
    influenceMatrix :: [t] -> Matrix Number
    influencePotential :: t -> Number -> Number -> Number
    tangentialVelocities :: [t] -> Vector Number
    liftCoefficient :: [t] -> Number
    error :: [t] -> Number

consMapper (SourcePanel2D _ _) = Source2D
consMapper (DoubletPanel2D _ _) = Doublet2D
consMapper (VortexPanel2D _ _) = Vortex2D

integrand panel f a b x y s = f (consMapper panel a s b) x y

instance Fluid2D Solution where
    velocity (Source2D str x0 y0) x y = (str/(2*pi)*(x - x0)/((x - x0)^2 + (y - y0)^2), str/(2*pi)*(y - y0)/((x - x0)^2 + (y - y0)^2))
    velocity (Uniform2D mag ang) x y = (mag*cos ang, mag*sin ang)
    velocity (Doublet2D str x0 y0) x y = (-str/(2*pi)*((x - x0)^2 - (y - y0)^2)/((x - x0)^2 + (y - y0)^2)^2, -str/(2*pi)*2*(x - x0)*(y - y0)/((x - x0)^2 + (y - y0)^2)^2)
    velocity (Vortex2D str x0 y0) x y = (-str/(2*pi)*(y - y0)/((x - x0)^2 + (y - y0)^2), str/(2*pi)*(x - x0)/((x - x0)^2 + (y - y0)^2))

    potential (Source2D str x0 y0) x y = str/(4*pi)*log ((x - x0)^2 + (y - y0)^2)
    potential (Uniform2D mag ang) x y = mag*(x*cos ang + y*sin ang)
    potential (Doublet2D str x0 y0) x y = -str/(2*pi)*(x - x0)/((x - x0)^2 + (y - y0)^2)
    potential (Vortex2D str x0 y0) x y = str/(2*pi)*atan2 (y - y0) (x - x0)

    stream (Source2D str x0 y0) x y = str/(2*pi)*atan2 (y - y0) (x - x0)
    stream (Uniform2D mag ang) x y = mag*(y*cos ang - x*sin ang)
    stream (Doublet2D str x0 y0) x y = -str/(2*pi)*(y - y0)/((x - x0)^2 + (y - y0)^2)
    stream (Vortex2D str x0 y0) x y = -str/(4*pi)*log ((x - x0)^2 + (y - y0)^2)

-- instance Fluid2D Panel where
    

instance SolvePanel Panel where
    influencePotential panel x y = (integral int 0 len) where
        int = integrand panel potential 1 0 x' y'
        (x', y') = rotation x y angle
        len = panelLength panel
        angle = panelAngle panel

    influenceVelocity panel x y = invRotation u' v' angle where
        u' = integral (fst . int) 0 len
        v' = integral (snd . int) 0 len
        int = integrand panel velocity 1 0 x' y'
        (x', y') = rotation x y angle
        len = panelLength panel
        angle = panelAngle panel

    influenceMatrix panels@(SourcePanel2D _ _ : panels') = a where 
        srcMat = 
            [ [ if panel_i == panel_j then 0.5*influencePotential panel_j xc yc else potential panel_j xc yc | panel_j <- panels, let (xc, yc) = colPoint (start panel_i) (end panel_i) ] | panel_i <- panels ]

    influenceMatrix panels@(DoubletPanel2D _ _ : panels') bc = 
        | map toLower bc == "dirichlet" = fromLists [ [ if panel_i == panel_j then 0.5 else influencePotential panel_j xc yc | panel_j <- panels, let (xc, yc) = colPoint (start panel_i) (end panel_i) ] | panel_i <- panels ]
        | map toLower bc == "neumann" = a where
        wokeVector = [ influenceVelocity wokePanel xc yc, let (xc, yc) = colPoint (start panel) (end panel) | panel <- panels ]
        dubMat = [ [ if panel_i == panel_j then 0.5 else sum $ influenceVelocity panel_j xc yc | panel_j <- panels, let (xc, yc) = colPoint (start panel_i) (end panel_i) ] | panel_i <- panels ]
        kutta = 
        a = fromBlocks [ [ dubMat, col wokeVector ],
                         [ kutta ,       0        ] ]
                
    -- influenceMatrix panels@(VortexPanel2D _ _ : panels') = [ [ if panel_i == panel_j then 0.5*potential panel_j xc yc else potential panel_j xc yc | panel_j <- panels, let (xc, yc) = colPoint (start panel_i) (end panel_i) ] | panel_i <- panels ]

wokePanel panels = DoubletPanel2D (xs, ys) (100000*xs, ys) where 
    trailing = head panels
    (xs, ys) = start trailing

eNorm = sqrt . sum . map (^2)

colPoint (xs, ys) (xe, ye) = ( (xs + ys)/2, (xe + ye)/2 )

panelCoords x y xs ys angle = rotation (x - xs) (y - ys) angle
rotation x y angle = (x*cos angle + y*sin angle, -x*sin angle + y*cos angle)
invRotation x y angle = (x*cos angle - y*sin angle, x*sin angle + y*cos angle)

panelLength panel = eNorm [xe1 - xs1, ye1 - ys1] where 
    (xs1, xe1) = start panel
    (ys1, ye1) = end panel

panelAngle panel = atan2 (ye1 - ys1) (xe1 - xs1) where
    (xs1, xe1) = start panel
    (ys1, ye1) = end panel