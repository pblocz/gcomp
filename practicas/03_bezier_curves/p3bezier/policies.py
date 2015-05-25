#!/usr/bin/env python2
# coding=utf-8

'''
Computational Geometry Assignments | (c) 2015 Pablo Cabeza & Diego González
license: [modified BSD](http://opensource.org/licenses/BSD-3-Clause)
'''

import numpy as np
from scipy.special import binom


class BezierPolicy(object):
    '''
    Base class used to implement different curve
    methods. It computes the curve given by the polygon.
    '''

    def __init__(self, polygon):
        self.polygon = np.array(polygon)
        self.N = len(self.polygon)

    def __call__(self, npoints=100, *args, **kwargs):
        return self.compute(np.linspace(0, 1, npoints),
                            *args, **kwargs)


class BernsteinPolicy(BezierPolicy):
    '''Naive Bernstain computation using numpy'''
    def _berstein(self, t):
        res = np.zeros((self.N + 1, len(t),))
        for i in range(self.N):
            res[i, :] = binom(self. N - 1, i) * t**i * (1 - t)**(self.N - i - 1)
        return res

    def compute(self, t):
        bernstein = self._berstein(t)
        curve_x = sum(self.polygon[i, 0] * bernstein[i, :]
                      for i in range(self.N))
        curve_y = sum(self.polygon[i, 1] * bernstein[i, :]
                      for i in range(self.N))
        return np.vstack((curve_x, curve_y,)).T


class FastBernsteinPolicy(BezierPolicy):
    '''Optimized Berstein computation (more vectorized), that uses precomputed binom'''

    from p3bezier.data import genen_binom as binom
    table = binom.table

    @staticmethod
    def _fastexp(t, n):
        '''
        Create a table of (n+1)xlen(t) elements where:

            table[i] = t**i

        It makes extensive use of numpy views to improve performance.
        It does O(log N) calls to numpy
        '''

        r = np.empty([n + 1, len(t)])
        r[0, :].fill(1); r[1, :] = t

        pi, i = 1, 1
        while i < n:
            newi = min(2 * i + 1, n)
            np.multiply(r[i, :], r[pi:pi + (newi - i)],  # multiply this two
                        r[i + 1:newi + 1, :])  # output array
            pi = i; i = newi
        return r

    @classmethod
    def _binom(cls, n, k):
        '''
        binomial implementation using a table cache, preinitialised
        '''
        if n == -1 or k == -1: return 0
        if k > n - k: k = n - k
        index = ((n * (n + 1)) >> 1) + k

        if len(cls.table) < index + k:
            cls.table.extend([-1] * ((((n + 1) * (n + 2)) >> 1) - len(cls.table)))
        if cls.table[index] == -1:
            cls.table[index] = cls._binom(n - 1, k) + cls._binom(n - 1, k - 1)

        return cls.table[index]

    @classmethod
    def _binom_row(cls, n):
        '''Returns the whole row on binomial coefficients of n'''

        base_idx = ((n * (n + 1)) >> 1)
        if len(cls.table) < base_idx + n + 1:
            return np.array(cls._binom(n, i) for i in range(base_idx, base_idex + n + 1))

        return np.array(cls.table[base_idx:base_idx + n + 1])

    def _berstein(self, t):
        e = self._fastexp(t, self.N - 1)
        res = e * e[:: - 1, :: - 1]
        np.multiply(res.T, self._binom_row(self.N - 1), res.T)
        return res

    def compute(self, t):
        bernstein = self._berstein(t)
        return np.dot(bernstein.T, self.polygon)


class FasterBernsteinPolicy(FastBernsteinPolicy):
    '''It uses pre-computed coefficients, currently there are up to degree 30'''
    from p3bezier.data.berncoef import coefficients_np as coefficients

    def compute(self, t):
        e = self._fastexp(t, self.N - 1)
        coef = np.dot(self.coefficients[self.N - 2], self.polygon)
        return np.dot(coef.T, e).T


class DeCasteljauPolicy(BezierPolicy):
    '''Compute Bernstein polynomial evaluted in the given points using DeCasteljau algorithm'''

    def compute(self, t):
        prevx = np.dot(self.polygon[:, 0][None].T, t[None])
        actx = np.empty(prevx.shape)
        prevy = np.dot(self.polygon[:, 1][None].T, t[None])
        acty = np.empty(prevy.shape)
        for j in range(1, self.N):
            for i in range(self.N - j):
                actx[i] = prevx[i + 1] * t + prevx[i] * (1 - t)
                acty[i] = prevy[i + 1] * t + prevy[i] * (1 - t)

            prevx, actx = actx, prevx
            prevy, acty = acty, prevy

        return np.vstack((prevx[0, :], prevy[0, :])).T


class DeCasteljauFastPolicy(BezierPolicy):
    '''Faster DeCasteljau algorithm that uses more numpy'''

    def compute(self, t):
        rt = t[::-1]
        prevx = np.dot(self.polygon[:, 0][None].T, t[None])
        actx = np.empty(prevx.shape)
        prevy = np.dot(self.polygon[:, 1][None].T, t[None])
        acty = np.empty(prevx.shape)

        for j in range(1, self.N):
            np.multiply(np.roll(prevx, -1, axis=0), t, actx)
            np.add(actx, np.multiply(prevx, rt), actx)

            np.multiply(np.roll(prevy, -1, axis=0), t, acty)
            np.add(acty, np.multiply(prevy, rt), acty)

            prevx, actx = actx[:-1], prevx[:-1]
            prevy, acty = acty[:-1], prevy[:-1]

        return np.vstack((prevx[0, :], prevy[0, :])).T


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from p3cbezier import bernstein as cberns
    poly = [[0, 0], [1, 1.5], [2, 2], [1.25, 0.7]]

    ffbcurve = FasterBernsteinPolicy(poly)()
    fbcurve = FastBernsteinPolicy(poly)()
    bcurve = BernsteinPolicy(poly)()
    ccurve = np.array(cberns(poly, 100))

    dccurve = DeCasteljauPolicy(poly)()
    dfccurve = DeCasteljauFastPolicy(poly)()

    poly = np.array(poly)
    plt.figure(1)
    plt.subplot(421)
    plt.plot(bcurve[:, 0], bcurve[:, 1])
    plt.plot(poly[:, 0], poly[:, 1], "ro")

    plt.subplot(422)
    plt.plot(fbcurve[:, 0], fbcurve[:, 1])
    plt.plot(poly[:, 0], poly[:, 1], "ro")

    plt.subplot(423)
    plt.plot(ffbcurve[:, 0], ffbcurve[:, 1])
    plt.plot(poly[:, 0], poly[:, 1], "ro")

    plt.subplot(424)
    plt.plot(ccurve[:, 0], ccurve[:, 1])
    plt.plot(poly[:, 0], poly[:, 1], "ro")

    plt.subplot(425)
    plt.plot(dccurve[:, 0], dccurve[:, 1])
    plt.plot(poly[:, 0], poly[:, 1], "ro")

    plt.subplot(426)
    plt.plot(dfccurve[:, 0], dfccurve[:, 1])
    plt.plot(poly[:, 0], poly[:, 1], "ro")
    plt.show()

    def eval_bezier(degree, t, algo):
        P = np.random.uniform(-20, 20, (degree + 1, 2))
        curve = algo(P)
        return curve  # numpy array of size (num_points, 2)

    import timeit
    degree = 15
    num_points = 100
    number = 10000
    tt = np.linspace(0, 1, num_points)

    def berns(P): return BernsteinPolicy(P)(npoints = num_points)
    def fberns(P): return FastBernsteinPolicy(P)(npoints = num_points)
    def ffberns(P): return FasterBernsteinPolicy(P)(npoints = num_points)

    def cast(P): return DeCasteljauPolicy(P)(npoints = num_points)
    def fcast(P): return DeCasteljauFastPolicy(P)(npoints = num_points)

    def ccberns(P): return cberns(P, num_points)

    print(timeit.timeit("eval_bezier(degree, tt, fberns)",
                        setup="from __main__ import eval_bezier, tt, degree, fberns, FastBernsteinPolicy", number=number))

    print(timeit.timeit("eval_bezier(degree, tt, ffberns)",
                        setup="from __main__ import eval_bezier, tt, degree, ffberns, FasterBernsteinPolicy", number=number))

    print(timeit.timeit("eval_bezier(degree, tt, berns)",
                        setup="from __main__ import eval_bezier, tt, degree, berns, BernsteinPolicy", number=number))

    print(timeit.timeit("eval_bezier(degree, tt, ccberns)",
                        setup="from __main__ import eval_bezier, tt, degree, ccberns, FastBernsteinPolicy", number=number))

    print(timeit.timeit("eval_bezier(degree, tt, cast)",
                        setup="from __main__ import eval_bezier, tt, degree, cast, DeCasteljauPolicy", number=number))

    print(timeit.timeit("eval_bezier(degree, tt, fcast)",
                        setup="from __main__ import eval_bezier, tt, degree, fcast, DeCasteljauFastPolicy", number=number))
