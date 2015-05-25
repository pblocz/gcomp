#!/usr/bin/env python2
# -*- encoding: utf-8 -*-

'''
Computational Geometry Assignments | (c) 2015 Pablo Cabeza & Diego González
license: modified BSD (http://opensource.org/licenses/BSD-3-Clause)
'''

import numpy as np


def quickhull(pts, segment=None):
    '''
    Computes the convex hull of a set of points.
    Inpired in [this quickhull](http://en.literateprograms.org/Quickhull_(Python,_arrays))

    To compute `dists`, we find a normal vector to the segment and then
    project each point to that vector using `numpy.dot` product.

    Arguments
    ---------
    - `pts`: `numpy.array(shape=(n,2))` of points to compute the hull

    Keyword argumens
    ----------------
    - `segment`: a pair (a,b) that divides the set of points. If `None`, then
      use as segment the extremes x-value points and initialize quickhull

    Return
    ------
    A list of point pairs that conforms the convex hull in clockwise order
    starting from the lowest lexicographic-wise point.
    '''

    if segment is None:  # initialization, find extremes and call quickhull
        pts = np.array(pts)
        xaxis = pts[:, 0]
        # m, M = pts[np.argmin(xaxis)], pts[np.argmax(xaxis)]

        m = pts.repeat(xaxis == np.amin(xaxis), axis=0)
        m = m[np.argmin(m[:, 1])]

        M = pts.repeat(xaxis == np.amax(xaxis), axis=0)
        M = M[np.argmax(M[:, 1])]

        return quickhull(pts, (m, M)) + quickhull(pts, (M, m))[1:]

    sB, sE = segment; vector = sE - sB
    dists = np.dot(pts - sB, [-vector[1], vector[0]])
    outer = pts.repeat(dists > 0, axis=0)  # take dists > eps elements
    # outer = pts[dists > eps]

    if not outer.size: return segment  # base case: no points above

    # recursive: find pivot and repeat over new segments
    pivot = pts[np.argmax(dists)]
    return quickhull(outer, segment=(sB, pivot)) +\
           quickhull(outer, segment=(pivot, sE))[1:]


def quickhull_recursion(pts, segment):
    '''
    Given a set of points and a segment dividing them, compute the *upper*
    convex hull.

    See `quickhull.quickhull` for more details.
    '''
    sB, sE = segment; vector = sE - sB
    dists = np.dot(pts - sB, [-vector[1], vector[0]])
    outer = pts.repeat(dists > 0, axis=0)  # take dists > eps elements
    # outer = pts[dists > eps]

    if not outer.size: return segment  # base case: no points above

    # recursive: find pivot and repeat over new segments
    pivot = pts[np.argmax(dists)]
    return quickhull_recursion(outer, segment=(sB, pivot)) +\
           quickhull_recursion(outer, segment=(pivot, sE))[1:]


def convex_hull(pts):
    '''
    Compute the convex hull of a set of points using quickhull algorithm.
    It initializes the algorithm and then calls `quickhull.quickhull_recursion`.

    See `quickhull.quickhull` for a more complete, *friendly*, version. It is
    split in 2 functions for performance.
    '''
    pts = np.array(pts)
    xaxis = pts[:, 0]

    m = pts.repeat(xaxis == np.amin(xaxis), axis=0)
    m = m[np.argmin(m[:, 1])]

    M = pts.repeat(xaxis == np.amax(xaxis), axis=0)
    M = M[np.argmax(M[:, 1])]

    return quickhull_recursion(pts, (m, M)) +\
           quickhull_recursion(pts, (M, m))[1:]
