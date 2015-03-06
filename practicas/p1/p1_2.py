# -*- coding: utf-8 -*-
'''
Practica 0

@author: Pablo Cabeza García
'''
import sympy as sp, numpy as np
from sympy import * #solve, symbols, Symbol, cos, sqrt, sin, lambdify, solve, pi, cosh, sinh
from sympy.interactive import printing
from sympy.parsing.sympy_parser import parse_expr
from sympy.plotting import plot_parametric as pl_para
from scipy.spatial import distance
from numpy import arange
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

printing.init_printing(pretty_print=True, use_latex=True)

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


def signature(Y):
    x,y = Y
    dx,dy = x.diff(t), y.diff(t)
    d2x,d2y = x.diff(t,2), y.diff(t,2)     


    norma_dY = sp.sqrt(dx**2 + dy**2)
    k = (dx * d2y - dy * d2x) / norma_dY**3 ; k=k.simplify()
    k_s = k.diff(t) / norma_dY; k_s = k_s.simplify()

    print (k,k_s)

    return (k,k_s)
    

def signature_check(Y1, I1, Y2, I2, eps = 0.1, step=0.01,plot=False):
    k1,ks1 = [lambdify(s.free_symbols or [t],s,"numpy") for s in signature(Y1)]
    k2,ks2 = [lambdify(s.free_symbols or [t],s,"numpy") for s in signature(Y2)]

    r1 = arange(I1[0],I1[1], step); r2 = arange(I2[0],I2[1], step)

    print type(r1)
    # print k1(r1)
    # print ks1(r1)

    def transel(el,r):
        if isinstance(el,np.ndarray): return el
        elif not isinstance(el,list): return [el]*len(r)
        return el
        
    ak1 = [transel(el,r1) for el in [k1(r1),ks1(r1)]]
    ak2 = [transel(el,r2) for el in [k2(r2),ks2(r2)]]

    p1 = zip(*ak1)
    p2 = zip(*ak2)

    dist = distance.cdist(p1,p2)
    dist_b = [[j < eps for j in i] for i in dist]
    trues = sum((i.count(True)) for i in dist_b)
    length = len(p1)*len(p2)
    
    print trues, length, 100*float(trues)/length

    if plot:
        # pl_para(ak1[0],ak1[1],(t,I1[0],I1[1]))
        # pl_para(ak2[0],ak2[1],(t,I2[0],I2[1]))

        plt.axes().set_aspect('equal')

        curve = plt.plot(ak1[0],ak1[1])
        add_arrow_to_line2D(plt.axes(), curve)
        
        plt.show()
        curve = plt.plot(ak2[0],ak2[1])
        add_arrow_to_line2D(plt.axes(), curve)
        plt.show()



    # nb=1
    # lenj=len(dist_b[0])
    # for i in range(len(dist_b)):
    #     l = []
    #     for j in range(lenj):
    #         i0,j0 = i, j
    #         def cbanda(i,j,nb=nb): 
    #             banda = dist_b[i][max(0,j-nb):min(j+nb+1,lenj)]
    #             return banda.count(True) != 0
                
    #         while cbanda(i0,j0): 
    #             i0+=1; j0+=1
    #         if (i,j) != (i0,j0): print (i,j), i0-i

if __name__ == '__main__':

    t, x, y = sp.var('t, x, y')

    # signature_check((t-1,t), t , (0,1), 
    #                 (2*t -2,3-t),t,(-1,0))

    # i1 = (t - 1,t); i1I=(0,1)
    # i2 = (2*t - 5, 3 - t); i2I = (-1,0)
    # signature_check(i1,i1I, i2,i2I)

    # ii1 = (2*cos(t),3*sin(t)); ii1I = (0,2*np.pi)
    # ii2 = (3*cos(t) ,2*sin(t)); ii2I = (0,2*np.pi)
    # signature_check(ii1,ii1I, ii2,ii2I)

    # iii1 = (t, 1  /(2*t)); iii1I=(0.1,10)
    # iii2 = (cosh(t),sinh(t)); iii2I = (0,1)
    # signature_check(iii1,iii1I, iii2,iii2I)

    # iv1 = (t,t**2); iv1I = (-2,2)
    # iv2 = ((-1/2.0)*sqrt(3)*log(t)**2 + (1/2.0)*log(t) + 1,
    #        ((1/2.0)*log(t)**2 + (1/2.0)*sqrt(3)*log(t) - 1))
    # iv2I = (0.1,10)
    # signature_check(iv1,iv1I, iv2,iv2I)

    i1 = (t,sin(t)); i1I=(10,40)
    i2 = (50-t,sin(50-t)); i2I=(10,40)
    signature_check(i1,i1I, i2,i2I,plot=True,step=0.1)
    
    
    or(t)
    or'(t)

    or(N-t)
    or'(N-t)*(-1)


    (dx(t) * d2y(t) - dy(t) * d2x(t))
