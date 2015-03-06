from scipy.integrate import odeint # instalarla
import matplotlib.pyplot as plt
import numpy as np
from math import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

def add_arrow_to_line2D(
    axes, line, arrow_locs=[0.2, 0.4, 0.6, 0.8],
    arrowstyle='-|>', arrowsize=1, transform=None):
    """
    Add arrows to a matplotlib.lines.Line2D at selected locations.

    Parameters:
    -----------
    axes: 
    line: list of 1 Line2D obbject as returned by plot command
    arrow_locs: list of locations where to insert arrows, % of total length
    arrowstyle: style of the arrow
    arrowsize: size of the arrow
    transform: a matplotlib transform instance, default to data coordinates

    Returns:
    --------
    arrows: list of arrows
    """
    if (not(isinstance(line, list)) or not(isinstance(line[0], 
                                           mlines.Line2D))):
        raise ValueError("expected a matplotlib.lines.Line2D object")
    x, y = line[0].get_xdata(), line[0].get_ydata()

    arrow_kw = dict(arrowstyle=arrowstyle, mutation_scale=10 * arrowsize)

    color = line[0].get_color()
    use_multicolor_lines = isinstance(color, np.ndarray)
    if use_multicolor_lines:
        raise NotImplementedError("multicolor lines not supported")
    else:
        arrow_kw['color'] = color

    linewidth = line[0].get_linewidth()
    if isinstance(linewidth, np.ndarray):
        raise NotImplementedError("multiwidth lines not supported")
    else:
        arrow_kw['linewidth'] = linewidth

    if transform is None:
        transform = axes.transData

    arrows = []
    for loc in arrow_locs:
        s = np.cumsum(np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2))
        n = np.searchsorted(s, s[-1] * loc)
        arrow_tail = (x[n], y[n])
        arrow_head = (np.mean(x[n:n + 2]), np.mean(y[n:n + 2]))
        p = mpatches.FancyArrowPatch(
            arrow_tail, arrow_head, transform=transform,
            **arrow_kw)
        axes.add_patch(p)
        arrows.append(p)
    return arrows

def solve(k,I,s0,delta=0.1,init_cond=(0,0,1,0)):
    '''
    @arg k: funcion de curvatura
    @arg I: intervalo [x0,xN]
    @arg s0:  punto donde se dan las cond iniciales

    @kwarg delta: para el intervalo
    @kwarg init_cond: [x,y,dx,dy] iniciales
    '''

    def rhs_eqs(Y,s):
        x, y , dx, dy = Y
        return [dx, dy, -k(s) * dy, k(s) * dx] 

    I0,In = I
    intervs = [np.arange(s0,I0-delta,-delta), np.arange(s0,In + delta, delta)]
    solvs = [odeint(rhs_eqs, init_cond, interv) for interv in intervs]
    solv = np.concatenate((solvs[0][::-1],solvs[1]),axis=0)
    

    print intervs
    print intervs[0]
    # print [k(s) for s in intervs[0]]
    
    curve_x = solv[:,0]; curve_y = solv[:,1]

    plt.axes().set_aspect('equal')
    curve = plt.plot(curve_x, curve_y)
    add_arrow_to_line2D(plt.axes(), curve)
    plt.show()

if __name__ == "__main__":
    # solve(lambda s: 1, [0, 2*np.pi], np.pi)
    # solve(lambda s: s * (s - 3) * (s + 5), [0,10], 0, delta=0.01)
    solve(lambda x: 1/x, [10000,30000], 10000,init_cond=[0,0,0.5,0])
    #solve(lambda s: 4*cos(s)**2 + sin(s)**2, [-100,100], 0)
