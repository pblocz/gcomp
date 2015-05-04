#! /usr/bin/env python2
# -*- encoding: utf-8 -*-

'''
Different implementations for interpolating polynomials.

authors: Pablo Cabeza & Diego González
date: 29/04/2015
'''
import numpy as np
import matplotlib.pyplot as plt


def boor_coef(i, r, k, knots, t):
    denom = float(knots[i + k + 1 - r] - knots[i])
    if denom == 0: return 0
    return (t - knots[i]) / denom


def boor_step(i, r, k, knots, points, t):
    if r == 0:
        ret = np.empty_like(t); ret.fill(points[i+1])
        return ret

    coef = boor_coef(i, r, k, knots, t)
    boor1 = boor_step(i - 1, r - 1, k, knots, points, t)
    boor2 = boor_step(i, r - 1, k, knots, points, t)
    res = (1 - coef) * boor1 +\
            coef * boor2
    return res



def spline2d(a, b, psi, k, nu, A, num_dots):
    '''
    Computes a plane spline curve of order k defined on the interval
    [a, b] with knots psi, multiplicities nu and coefficients A.

    Parameters:
      - a, b -- ends of the interval, real numbers
      - psi -- list of breakpoints, a < psi[0] < ... < psi[-1] < b
      - k -- order of the curve, the degree is <= k - 1
      - nu -- list of integer multiplicities of each breakpoint,
        len(psi) = len(nu), 1 <= nu[i] < k
      - A -- list of coefficients of the B-spline basis,
        A = [[x0, y0], [x1, y1],..., [xN, yN]
      - num_dots -- number of dots of the spline to be plotted,
        uniformly spaced along the interval [a, b]

    Returns:
    the spline curve as a numpy array of size (2, num_dots)
    '''

    psi = np.array(psi, dtype = 'f')
    nu = np.array(nu, dtype = 'f')
    A = np.array(A, dtype = 'f')

    if len(nu)!=len(psi): return -1

    t = np.linspace(a, b, num_dots)

    nknots = 2*k + k*len(nu) - sum(nu)
    knots, k_points = np.empty(nknots), np.empty([nknots, 2])

    knots[:k].fill(a); knots[-k:].fill(b)
    k_points[:k][:, 0].fill(A[0, 0]); k_points[:k][:, 1].fill(A[0, 1])
    k_points[-k:][:, 0].fill(A[-1, 0]); k_points[-k:][:, 1].fill(A[-1, 1])

    idx = k
    for i, n in enumerate(nu):
        off = k - n
        knots[idx:idx + off].fill(psi[i])
        k_points[idx:idx + off][:, 0].fill(A[i, 0])
        k_points[idx:idx + off][:, 1].fill(A[i, 1])
        idx += off


    ret1, ret2, N = [], [], len(knots)
    for l in xrange(N):  # calculamos a_j^{k-1}
        if l == N - 1: break


        cond = (t >= knots[l]) & \
            ((t <= knots[l+1]) if l == N - k else (t < knots[l+1]))
        eval_points = t[cond]

        if len(eval_points) == 0: continue

        ret1.append(boor_step(l, k - 1, k - 1, knots, A[:,0], eval_points))
        ret2.append(boor_step(l, k - 1, k - 1, knots, A[:,1], eval_points))

    return np.array([np.concatenate(ret1), np.concatenate(ret2)])
