#!/usr/bin/env python2

'''
Computational Geometry Assignments | (c) 2015 Pablo Cabeza & Diego González
license: modified BSD
'''


import numpy as np
import matplotlib.pyplot as plt


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


# def andrews_hull2(pts):
#     pts = sorted(set(tuple(p) for p in pts))
#
#     def cross(O, A, B): return (A[0] - O[0]) * (B[1] - O[1]) - (A[1] - O[1]) * (B[0] - O[0])
#
#     r = []; lwr = r
#     for p in pts:
#         while len(r) >= 2 and cross(r[-2], r[-1], p) <= 0: r.pop()
#         r.append(p)
#
#     r = []; upp = r
#     for p in reversed(pts):
#         while len(r) >= 2 and cross(r[-2], r[-1], p) <= 0: r.pop()
#         r.append(p)
#
#     return np.array((lwr[:-1] + upp)[::-1])
#
convex_hull = andrews_hull


'''
if __name__ == "__main__":
    # assert convex_hull([(i / 10, i % 10)
    #                 for i in range(100)]) == [(0, 0), (9, 0), (9, 9), (0, 9)]

    npoints = 30

    # points = np.random.randint(-10, 10, size=(npoints, 2))
    # hull = np.array(quickhull(points))

    # points = np.array([[0, 0], [1, 0]])
    # points = np.array([[-1, 0], [1, 1], [1, 0], [-1, 0], [1, 0], [1, 1]])
    # points = np.array([[0, 0], [0, 1], [0, 2], [0, 3], [0, 1]])

    execfile("data2.py")
    points = np.array(points)

    hull = (convex_hull(points))
    print(hull)
    print(result)

    hull = np.array(hull)
    plt.plot(hull[:, 0], hull[:, 1])
    plt.plot(points[:, 0], points[:, 1], 'r.')
    plt.show()
'''
