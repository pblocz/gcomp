#!/usr/bin/env python2
# coding=utf-8

'''
Computational Geometry Assignments | (c) 2015 Pablo Cabeza & Diego González
license: [modified BSD](http://opensource.org/licenses/BSD-3-Clause)
'''

import sympy as sp
import numpy as np
from sympy import *  # solve, symbols, Symbol, cos, sqrt, sin, lambdify, pi...
from sympy.interactive import printing
from sympy.parsing.sympy_parser import parse_expr
from sympy.plotting import plot_parametric as pl_para
from scipy.spatial import distance
from numpy import arange
import matplotlib.pyplot as plt


def signature(Y):
    '''
    Compute the signature of a 2D curve.

    - `Y` = (x(t), y(t)); `sympy` expressions depending of a global symbol `t`

    **Return**: a pair (curvature, curvature_derivate) of `Y`
    '''
    x, y = Y
    dx, dy = x.diff(t), y.diff(t)
    d2x, d2y = x.diff(t, 2), y.diff(t, 2)

    norma_dY = sp.sqrt(dx**2 + dy**2)
    k = (dx * d2y - dy * d2x) / norma_dY**3; k = k.simplify()
    k_s = k.diff(t) / norma_dY; k_s = k_s.simplify()

    return (k, k_s)


def signature_check(Y1, I1, Y2, I2,
                    eps=0.1, step=0.01, ths=0.1, plot=False):
    '''
    Comprueba si la signatura de Y1 en el intervalo I1;
    Y2 en el intervalo I2 son iguales.

    kwargs
    ------
    - eps: epsilon con el que comprobar los ceros
    - step: el paso usado para calcular los puntos de los intervalos I1,I2
    - ths: el umbral a partir del cual considerar las curvas iguales
        (el porcentaje de puntos iguales en la matriz de distancias)
    - plot: si es True, entonces dibuja las signaturas
    '''
    k1, ks1 = [lambdify(s.free_symbols or [t], s, "numpy")
               for s in signature(Y1)]
    k2, ks2 = [lambdify(s.free_symbols or [t], s, "numpy")
               for s in signature(Y2)]

    r1 = arange(I1[0], I1[1], step)
    r2 = arange(I2[0], I2[1], step)

    # por alguna razón al aplicar el resultado de lambdify a una función que da
    # un mismo resultado siempre (por ejemplo idénticamente 0), en lugar de un
    # array devuelve un solo entero... Hemos hecho esto para resolverlo
    def transel(el, r):
        if isinstance(el, np.ndarray): return el
        elif not isinstance(el, list): return [el] * len(r)
        return el

    # Aplicamos la solución a las signaturas
    ak1 = [transel(el, r1) for el in [k1(r1), ks1(r1)]]
    ak2 = [transel(el, r2) for el in [k2(r2), ks2(r2)]]

    # Generamos los pares [(k,k_s)]
    p1 = zip(*ak1); p2 = zip(*ak2)

    # Calculamos las distancias y el número de elementos que están cerca
    dist = distance.cdist(p1, p2)
    dist_b = [[j < eps for j in i] for i in dist]
    trues = sum((i.count(True)) for i in dist_b)
    length = len(p1) * len(p2)

    # Para ayudar a visualizarlo añadimos el parámetro plot
    if plot:
        plt.plot(ak1[0], ak1[1])
        plt.show()
        plt.plot(ak2[0], ak2[1])
        plt.show()

    print("[%s] y [%s] iguales? %s" % (Y1, Y2, float(trues) / length >= ths))
    return float(trues) / length >= ths


if __name__ == '__main__':

    t, x, y = sp.var('t, x, y')

    iii1 = (t, 1 / (2 * t))
    iii1I = (0.1, 10)
    iii2 = (cosh(t), sinh(t))
    iii2I = (0, 1)
    signature_check(iii1, iii1I, iii2, iii2I)

    iv1 = (t, t**2)
    iv1I = (-2, 2)
    iv2 = ((-1 / 2.0) * sqrt(3) * log(t)**2 + (1 / 2.0) * log(t) + 1,
           ((1 / 2.0) * log(t)**2 + (1 / 2.0) * sqrt(3) * log(t) - 1))
    iv2I = (0.1, 10)
    signature_check(iv1, iv1I, iv2, iv2I)

    i1 = (exp(2 * t), 2 * exp(t))
    i1I = (-10, 10)
    signature_check(i1, i1I, i1, i1I, plot=True, step=0.1)
