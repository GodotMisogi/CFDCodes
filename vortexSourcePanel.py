import numpy as np
import sympy
from scipy import integrate
import time, sys, os
from matplotlib import pyplot as plt
from matplotlib import rc, cm
import math
from mpl_toolkits.mplot3d import Axes3D
import timeit

import importlib
import aeroClasses

importlib.reload(aeroClasses)
from aeroClasses import *

rc('font',**{'family':'serif'})
rc('text', usetex=True)

naca0012_filepath = os.path.join('resources', 'ClarkY.dat')
with open (naca0012_filepath, 'r') as file_name:
    x_coords, y_coords = np.loadtxt(file_name, dtype=float, delimiter='    ', unpack=True)

grid_size = 50
x_domain, y_domain = (-1, 2), (-0.75, 0.75)
x_dom, y_dom = np.linspace(x_domain[0], x_domain[1], grid_size), np.linspace(y_domain[0], y_domain[1], grid_size)
X, Y = np.meshgrid(x_dom, y_dom)

uniform_mag, uniform_ang = 1.0, 5.0
N_panels = 40
panels = cosine_panels(x_coords, y_coords, N_panels)
source_panels = [ SourcePanel2D(panel.xs, panel.ys, panel.xe, panel.ye) for panel in panels ]
uniform = Uniform2D(uniform_mag, uniform_ang)

vortex_naca0012 = VortexSourcePanel2DSolver(source_panels, uniform)
strengths = vortex_naca0012.solve_strengths()
vts = vortex_naca0012.tangential_velocity()
vels = vortex_naca0012.velocity(X, Y)
pots = vortex_naca0012.potential(X, Y)
cp = pressure_coefficient2D(vels[0], vels[1], uniform_mag)
error = vortex_naca0012.error()
cl = vortex_naca0012.lift_coefficient()

width = 10
fig22 = plt.figure(1, dpi=300)
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.plot(x_coords, y_coords, linestyle='--', linewidth=2)
plt.plot(np.append([panel.xs for panel in panels], panels[0].xs),
            np.append([panel.ys for panel in panels], panels[0].ys),
            linestyle='-', linewidth=1, marker='.', markersize=6, color='#CD2305')
plt.axis('scaled', adjustable='box')
plt.xlim(-0.1, 1.1)
plt.ylim(min([panel.ys for panel in panels]) - 0.05, max([panel.ys for panel in panels]) + 0.05)
plt.grid()

width = 8.0
height = (y_domain[1] - y_domain[0]) / (x_domain[1] - x_domain[0]) * width
fig23 = plt.figure(2, figsize=(width, height), dpi=300)
plt.grid()
plt.xlabel('$x$')
plt.ylabel('$C_{p_{local}}$')
plt.plot([panel.xc for panel in source_panels if panel.loc == 'upper'],
            [panel.cp for panel in source_panels if panel.loc == 'upper'],
            label='upper surface',
            color='r', linestyle='-', marker='.')
plt.plot([panel.xc for panel in source_panels if panel.loc == 'lower'],
            [panel.cp for panel in source_panels if panel.loc == 'lower'],
            label= 'lower surface',
            color='b', linestyle='-', marker='.')
plt.ylim(max([p.cp for p in source_panels]) + 0.05, min([p.cp for p in source_panels]) - 0.05)
plt.legend()

width = 8.0
height = (y_domain[1] - y_domain[0]) / (x_domain[1] - x_domain[0]) * width
fig24 = plt.figure(3, figsize=(width*1.2, height), dpi=300)
plt.streamplot(X, Y, vels[0], vels[1], density=3, zorder=2)
plt.contourf(X, Y, cp, 100, cmap=cm.rainbow)
plt.colorbar(label=r'$C_p$')
plt.plot(x_coords, y_coords,
            label='Coordinates',
            color='k', linestyle='--', linewidth=1)
plt.plot([p.xs for p in source_panels], [p.ys for p in source_panels], label='Panels')
plt.fill([p.xs for p in source_panels], [p.ys for p in source_panels], color='k', zorder=3)
plt.title(f'Number of vortex-source panels: {N_panels}')
plt.xlim((x_domain))
plt.ylim(y_domain);

plt.show()