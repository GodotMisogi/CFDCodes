import Math.GaussianQuadratureIntegration
import Numeric.LinearAlgebra

integral = nIntegrate256

type Number = Float

class Fluid2D t where
    velocity :: t -> Number -> Number -> (Number, Number)
    potential :: t -> Number -> Number -> Number
    stream :: t -> Number -> Number -> Number

-- class SolvePanel where
--     solveStrengths :: t1 -> Vector Number
--     tangentialVelocities :: t1 -> Vector Number
--     liftCoefficient :: t1 -> Number
--     error :: t1 -> Number
    
data Solution = 
    Source2D Number Number Number |
    Uniform2D Number Number |
    Doublet2D Number Number Number |
    Vortex2D Number Number Number 
    deriving (Eq, Show)

instance Fluid2D Solution where
    velocity (Source2D str x0 y0) x y = (str/(2*pi)*(x - x0)/((x - x0)^2 + (y - y0)^2), str/(2*pi)*(y - y0)/((x - x0)^2 + (y - y0)^2))
    velocity (Uniform2D mag ang) x y = (mag*cos ang, mag*sin ang)
    velocity (Doublet2D str x0 y0) x y = (-str/(2*pi)*((x - x0)^2 - (y - y0)^2)/((x - x0)^2 + (y - y0)^2)^2, -str/(2*pi)*2*(x - x0)*(y - y0)/((x - x0)^2 + (y - y0)^2)^2)
    velocity (Vortex2D str x0 y0) x y = (-str/(2*pi)*(y - y0)/((x - x0)^2 + (y - y0)^2), str/(2*pi)*(x - x0)/((x - x0)^2 + (y - y0)^2))

    potential (Source2D str x0 y0) x y = str/(4*pi)*log ((x - x0)^2 + (y - y0)^2)
    potential (Uniform2D mag ang) x y = mag*(x*cos ang + y*sin ang)
    potential (Doublet2D str x0 y0) x y = -str/(2*pi)*(x - x0)/((x - x0)^2 + (y - y0)^2)
    potential (Vortex2D str x0 y0) x y = str/(2*pi)*atan (y - y0)/(x - x0)

    stream (Source2D str x0 y0) x y = str/(2*pi)*atan (y - y0)/(x - x0)
    stream (Uniform2D mag ang) x y = mag*(y*cos ang - x*sin ang)
    stream (Doublet2D str x0 y0) x y = -str/(2*pi)*(y - y0)/((x - x0)^2 + (y - y0)^2)
    stream (Vortex2D str x0 y0) x y = -str/(4*pi)*log ((x - x0)^2 + (y - y0)^2)

data Panel =
    SourcePanel2D { start :: (Number, Number), end :: (Number, Number) } |
    DoubletPanel2D { start :: (Number, Number), end :: (Number, Number) } |
    VortexPanel2D { start :: (Number, Number), end :: (Number, Number) }
    deriving (Eq, Show)

instance Fluid2D Panel where
    -- velocity panel (x, y) = invRotation (u', v')  where
    --     u' = integral u 0 length $ rotation (x, y) angle
    --     v' = integral v 0 length $ rotation (x, y) angle
    --     (u, v) = velocity (Source2D 1 x0 0) x
    --     length = panelLength panel
    --     angle = panelAngle panel
    --     panel = SourcePanel2D start end

    -- velocity panel (x, y) = invRotation (u', v')  where
    --     u' = integral u 0 length $ rotation (x, y) angle
    --     v' = integral v 0 length $ rotation (x, y) angle
    --     (u, v) = velocity (Doublet2D 1 x0 0) x
    --     length = panelLength panel
    --     angle = panelAngle panel
    --     panel = DoubletPanel2D start end

    -- velocity panel (x, y) = invRotation (u', v')  where
    --     u' = integral u 0 length $ rotation (x, y) angle
    --     v' = integral v 0 length $ rotation (x, y) angle
    --     (u, v) = velocity (Vortex2D 1 x0 0) x
    --     length = panelLength panel
    --     angle = panelAngle panel
    --     panel = VortexPanel2D start end

    potential panel@(SourcePanel2D start end) x y = (integral integrand 0 len) where
        (x', y') = rotation (x, y) angle
        integrand = potential (Source2D 1 0 0) x
        len = panelLength panel
        angle = panelAngle panel

    -- potential panel@(DoubletPanel2D start end) (x, y) = (integral integrand 0 len) x' y' where
    --     (x', y') = rotation (x, y) angle
    --     integrand = potential (Doublet2D 1 0 0) x
    --     len = panelLength panel
    --     angle = panelAngle panel


    -- potential panel@(VortexPanel2D start end) (x, y) = (integral integrand 0 len) x' y' where
    --     (x', y') = rotation (x, y) angle
    --     integrand = potential (Vortex2D 1 0 0) x
    --     len = panelLength panel
    --     angle = panelAngle panel

eNorm = sqrt . sum . map (^2)

rotation (x, y) angle = (x*cos angle + y*sin angle, -x*sin angle + y*cos angle)
invRotation (x, y) angle = (x*cos angle - y*sin angle, x*sin angle + y*cos angle)

panelLength panel = eNorm [xe1 - xs1, ye1 - ys1] where 
    (xs1, xe1) = start panel
    (ys1, ye1) = end panel

panelAngle panel = atan (ye1 - ys1)/(xe1 - xs1) where
    (xs1, xe1) = start panel
    (ys1, ye1) = end panel