#!/usr/bin/env python2
# -*- encoding: utf-8 -*-

'''
Computational Geometry Assignments | (c) 2015 Pablo Cabeza & Diego González
license: modified BSD (http://opensource.org/licenses/BSD-3-Clause)
'''

import numpy as np


def andrews_hull(pts):
    '''
    Compute the convex hull using andrews hull algorithm.

    The algorithm and some inspirational code can be found [here](http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain))

    Arguments
    ---------
    - `pts`: `numpy.array(shape=(n,2))` or list of points pairs to compute hull

    Return
    ------
    A list of point pairs that conforms the convex hull clockwise starting from
    the lowest lexicographic-wise point
    '''
    pts = sorted(set(tuple(p) for p in pts))

    # Faster than np.cross for single points, just work with lists
    def cross(O, A, B): return (A[0] - O[0]) * (B[1] - O[1]) - (A[1] - O[1]) * (B[0] - O[0])

    def dome(pts):
        r = []
        for p in pts:
            while len(r) >= 2 and cross(r[-2], r[-1], p) <= 0: r.pop()
            r.append(p)
        return r

    return (dome(pts)[:-1] + dome(reversed(pts)))[::-1]


convex_hull = andrews_hull
