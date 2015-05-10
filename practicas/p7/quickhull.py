#!/usr/bin/env python2

import numpy as np
import matplotlib.pyplot as plt


def pplot(points, *args, **kwargs):
    "plot points in the form of numpy.array(shape=(m,2))"
    plt.plot(points[:, 0], points[:, 1], *args, **kwargs)


def quickhull(pts, eps=0.01, segment=None):
    '''
    Computes the convex hull of a set of points.
    Inpired in [this quickhull](http://en.literateprograms.org/Quickhull_(Python,_arrays))

    To compute `dists`, we find a normal vector to the segment and then
    project each point to that vector using the `np.dot` product.

    Arguments
    ---------
    - `pts`: `numpy.array(shape=(n,2))` of points to compute the hull

    Keyword argumens
    ----------------
    - `eps`: the threadhold to compare with zero
    - `segment`: a pair (a,b) that divides the set of points. If `None`, then
      use as segment the extremes x-value points and initialize quickhull

    Return
    ------
    The resulting hull in counterclockwise order as a (m,2) shape `numpy.array`
    '''

    if segment is None:  # initialization, find extremes and call quickhull
        xaxis = pts[:, 0]

        m = pts.repeat(xaxis == np.amin(xaxis), axis=0)
        m = m[np.argmin(m[:, 1])]
        
        M = pts.repeat(xaxis == np.amax(xaxis), axis=0)
        M = M[np.argmax(M[:, 1])]

        return quickhull(pts, eps, (m, M)) + quickhull(pts, eps, (M, m))[1:]

    sB, sE = segment; vector = sE - sB
    dists = np.dot(pts - sB, [-vector[1], vector[0]])
    outer = pts.repeat(dists > eps, axis=0)  # take dists > eps elements

    if not outer.size: return segment  # base case: no points above

    # recursive: find pivot and repeat over new segments
    pivot = pts[np.argmax(dists)]
    return quickhull(outer, eps, segment=(sB, pivot)) +\
           quickhull(outer, eps, segment=(pivot, sE))[1:]


def convex_hull(pts): return quickhull(np.array(pts), eps=0.0)


'''
if __name__ == "__main__":

    npoints = 30

    # points = np.random.randint(-10, 10, size=(npoints, 2))
    # hull = np.array(quickhull(points))

    # points = np.array([[-1, 0], [1, 0], [1, 1], [-1, 0], [1, 0], [1, 1]])
    # points = np.array([[0, 0], [0, 1], [0, 2], [0, 3], [0, 1]])

    execfile("data.py")
    points = np.array(points)
    hull = np.array(convex_hull(points))
    print(hull)

    plt.plot(hull[:, 0], hull[:, 1])
    plt.plot(points[:, 0], points[:, 1], 'r.')

    # import time
    # for i in xrange(len(hull) - 1):
    #     pplot(np.vstack((hull[i], hull[i + 1], )))
    #     time.sleep(1)

    plt.show()
'''
