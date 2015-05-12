#!/usr/bin/env python2

import numpy as np
import matplotlib.pyplot as plt

# http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain

def convex_hull(points):
    """Computes the convex hull of a set of 2D points.

    Input: an iterable sequence of (x, y) pairs representing the points.
    Output: a list of vertices of the convex hull in counter-clockwise order,
      starting from the vertex with the lexicographically smallest coordinates.
    Implements Andrew's monotone chain algorithm. O(n log n) complexity.
    """
    points = [tuple(p) for p in points]

    # Sort the points lexicographically (tuples are compared lexicographically).
    # Remove duplicates to detect the case we have just one unique point.
    points = sorted(set(points))

    # Boring case: no points or a single point, possibly repeated multiple times.
    if len(points) <= 1:
        return points

    # 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
    # Returns a positive value, if OAB makes a counter-clockwise turn,
    # negative for clockwise turn, and zero if the points are collinear.
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    # Build lower hull
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    # Build upper hull
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    # Concatenation of the lower and upper hulls gives the convex hull.
    # Last point of each list is omitted because it is repeated at the beginning of the other list.
    return upper[::-1][:-1] + lower[::-1]  # upper[:-1]


# Example: convex hull of a 10-by-10 grid.

'''
if __name__ == "__main__":
    # assert convex_hull([(i / 10, i % 10)
    #                 for i in range(100)]) == [(0, 0), (9, 0), (9, 9), (0, 9)]

    npoints = 30

    # points = np.random.randint(-10, 10, size=(npoints, 2))
    # hull = np.array(quickhull(points))

    # points = np.array([[0, 0], [1, 0]])
    points = np.array([[-1, 0], [1, 0], [1, 1], [-1, 0], [1, 0], [1, 1]])
    # points = np.array([[0, 0], [0, 1], [0, 2], [0, 3], [0, 1]])

    # execfile("data.py")
    # points = np.array(points)

    hull = np.array(convex_hull(points))
    print(hull)

    plt.plot(hull[:, 0], hull[:, 1])
    plt.plot(points[:, 0], points[:, 1], 'r.')
'''