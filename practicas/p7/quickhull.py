#!/usr/bin/env python2

import numpy as np
import matplotlib.pyplot as plt


def pplot(points, *args, **kwargs):
    "plot points in the form of numpy.array(shape=(m,2))"
    plt.plot(points[:, 0], points[:, 1], *args, **kwargs)


def qhull(points):
    '''
    Idea sacada de (http://en.literateprograms.org/Quickhull_(Python,_arrays))
    '''
    def link(a, b): return np.concatenate((a, b[1:]))

    def edge(a, b): return np.concatenate(([a], [b]))

    def dome(points, base):
        h, t = base
        dists = np.dot(points - h, np.dot([(0, -1), (1, 0)], (t - h)))
        outer = np.repeat(points, dists > 0, axis=0)

        if len(outer):
            pivot = points[np.argmax(dists)]
            return link(dome(outer, edge(h, pivot)),
                        dome(outer, edge(pivot, t)))
        else:
            return base

    if len(points) > 2:
        axis = points[:, 0]
        base = np.take(points, [np.argmin(axis), np.argmax(axis)], axis=0)
        return link(dome(points, base),
                    dome(points, base[::-1]))
    else:
        return points


def quickhull(points, eps=0.01, segment=None):
    '''
    Computes the convex hull of a set of points.
    Inpired in [this quickhull](http://en.literateprograms.org/Quickhull_(Python,_arrays))

    Arguments
    ---------
    - `points`: `numpy.array(shape=(n,2))` of points to compute the hull

    Keyword argumens
    ----------------
    - `eps`: the threadhold to compare with zero
    - `segment`: a pair (m,M) that divides the set of points. If `None`, then
      use as segment the extremes x-value points and initialize quickhull

    Return
    ------
    The resulting hull in counterclockwise order as a (m,2) shape `numpy.array`
    '''

    # initialization
    if segment is None:
        xaxis = points[:, 0]
        m, M = points[np.argmin(xaxis)], points[np.argmax(xaxis)]
        return quickhull(points, eps, (m, M)) + quickhull(points, eps, (M, m))[1:]
    sB, sE = segment; vector = sE - sB

    dists = np.dot(points - sB, [-vector[1], vector[0]])
    outer = np.repeat(points, dists > eps, axis=0)  # faster than points[dists > 0] ?

    # base case: no points above
    if not outer.size: return segment

    # recursive: find pivot and repeat over new segments
    pivot = points[np.argmax(dists)]
    return quickhull(outer, eps, segment=(sB, pivot)) +\
           quickhull(outer, eps, segment=(pivot, sE))[1:]


if __name__ == "__main__":

    npoints = 30
    points = np.random.randint(-10, 10, size=(npoints, 2))

    hull = np.array(quickhull(points))
    print(hull)

    plt.plot(points[:, 0], points[:, 1], 'or')
    plt.plot(hull[:, 0], hull[:, 1])

    # import time
    # for i in xrange(len(hull) - 1):
    #     pplot(np.vstack((hull[i], hull[i + 1], )))
    #     time.sleep(1)

    plt.show()
