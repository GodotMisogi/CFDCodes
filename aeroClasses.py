import numpy as np

class Source2D:
    """
    Defines a source or sink.
    
    Parameters
    ---
    sigma: strength
    x_0: source x-location
    y_0: source y-location
    """
    
    def __init__(self, sigma, x_0, y_0):
        self.strength = sigma
        self.x_0 = x_0
        self.y_0 = y_0
        
    def velocity(self, x, y):
        u = self.strength/(2*np.pi)*(x - self.x_0)/((x - self.x_0)**2 + (y - self.y_0)**2)
        v = self.strength/(2*np.pi)*(y - self.y_0)/((x - self.x_0)**2 + (y - self.y_0)**2)
        
        return (u, v)
    
    def stream_func(self, x, y):
        psi = self.strength/(2*np.pi)*np.arctan2(y - self.y_0, x - self.x_0)
        
        return psi
        
    def potential(self, x, y):
        pot = self.strength/(4*np.pi)*np.log(np.sqrt((x - self.x_0)**2 + (y - self.y_0)**2))
        
        return pot
    
class Uniform2D:
    """
    Defines a uniform flow.
    
    Parameters:
    1. u_inf = magnitude
    2. alpha = angle
    """
    
    def __init__(self, u_inf, alpha):
        self.U = u_inf
        self.angle = alpha
        
    def velocity(self, x, y):
        u = self.U*(np.cos(np.deg2rad(self.angle)))
        v = self.U*(np.sin(np.deg2rad(self.angle)))
        
        return (u, v)
    
    def stream_func(self, x, y):
        psi = self.U*(y*np.cos(np.deg2rad(self.angle)) - x*np.sin(np.deg2rad(self.angle)))
        
        return psi
        
    def potential(self, x, y):
        pot = self.U*(x*np.cos(np.deg2rad(self.angle)) + y*np.sin(np.deg2rad(self.angle)))
        
        return pot

class Doublet2D:
    """
    Defines a doublet.
    
    Parameters
    ---
    kappa: strength
    x_0: doublet x-location
    y_0: doublet y-location
    """
    
    def __init__(self, kappa, x_0, y_0):
        self.strength = kappa
        self.x_0 = x_0
        self.y_0 = y_0
        
    def velocity(self, x, y):
        u = -self.strength/(2*np.pi)*((x - self.x_0)**2 - (y - self.y_0)**2)/((x - self.x_0)**2 + (y - self.y_0)**2)**2
        v = -self.strength/(2*np.pi)*2*(x - self.x_0)*(y - self.y_0)/((x - self.x_0)**2 + (y - self.y_0)**2)**2
        
        return (u, v)
    
    def stream_func(self, x, y):
        psi = -self.strength/(2*np.pi)*(y - self.y_0)/((x - self.x_0)**2 + (y - self.y_0)**2)
        
        return psi
        
    def potential(self, x, y):
        # CHECK THIS YOU RETARD
        pot = self.strength/(2*np.pi)*np.log(np.sqrt((x - self.x_0)**2 + (y - self.y_0)**2))
        
        return pot

class Vortex2D:
    """
    Defines a vortex.
    
    Parameters
    ---
    Gamma: strength
    x_0: vortex x-location
    y_0: vortex y-location
    """
    
    def __init__(self, gamma, x_0, y_0):
        self.strength = gamma
        self.x_0 = x_0
        self.y_0 = y_0
        
    def velocity(self, x, y):
        u = -self.strength/(2*np.pi)*(y - self.y_0)/((x - self.x_0)**2 + (y - self.y_0)**2)
        v = self.strength/(2*np.pi)*(x - self.x_0)/((x - self.x_0)**2 + (y - self.y_0)**2)
        
        return (u, v)
    
    def stream_func(self, x, y):
        psi = -self.strength/(4*np.pi)*np.log(np.sqrt((x - self.x_0)**2 + (y - self.y_0)**2))
        
        return psi
        
    def potential(self, x, y):
        pot = self.strength/(2*np.pi)*np.arctan2(y - self.y_0, x - self.x_0)
        
        return pot

class InfiniteVortices:
    """
    Defines an infinite row of vortices.
    
    Parameters
    ---
    Gamma: strength
    a: spacing
    """
    
    def __init__(self, gamma, a):
        self.strength = gamma
        self.a = a
        
    def velocity(self, x, y):
        u = -self.strength/(2*self.a)*np.sinh(2*np.pi*y/self.a)/(np.cosh(2*np.pi*y/self.a) - np.cos(2*np.pi*x/self.a))
        v = self.strength/(2*self.a)*np.sin(2*np.pi*x/self.a)/(np.cosh(2*np.pi*y/self.a) - np.cos(2*np.pi*x/self.a))
        
        return (u, v)


def joukowski(z, c):
    """
    Computes the Joukowski transformation of a complex number z.
    
    Parameters:
    1. z: complex number in Cartesian coordinates
    2. c: some parameter
    """
    xi = z + c**2/z
    return xi

def pressure_coefficient2D(u, v, u_inf):
    "Computes the pressure coefficient in 2 dimensions."
    cp = 1. - (u**2 + v**2)/u_inf**2
    return cp