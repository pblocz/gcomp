#!/usr/bin/env python
# coding=utf-8

'''
Computational Geometry Assignments | (c) 2015 Pablo Cabeza & Diego González
license: [modified BSD](http://opensource.org/licenses/BSD-3-Clause)
'''

# Standard python libraries
import sys
import logging
from collections import namedtuple

try: import Queue as Q  # ver. < 3.0
except ImportError: import queue as Q

# Custom Installed libraries
import matplotlib.pyplot as plt
from bintrees import RBTree as btree


logger = logging.getLogger(__name__)  # Get current app logger


class Point(namedtuple("Point_tuple", ['x', 'y'])):
    '''Pair of points ordered lexicophically with `+`, `-`, `*` overloaded'''

    def pad(o): return o if type(o) == Point else Point(*[o for i in range(2)])
    def __add__(s, o): return Point(*[a + b for a, b in zip(s, Point.pad(o))])
    def __sub__(s, o): return Point(*[a - b for a, b in zip(s, Point.pad(o))])

    def __mul__(s, o):
        '''if o is `int` then multiply each component, else dot product'''
        if isinstance(o, Point): return sum(a * b for a, b in zip(s, Point.pad(o)))
        else: return Point(*[a * b for a, b in zip(s, Point.pad(o))])

    def __repr__(s): return "Point(x=%.2f, y=%.2f)" % s


class Segment(namedtuple('Segment_tuple', ['left', 'right'])):
    '''A segment formed by a left and right `segment.Point`'''

    X = None
    '''Current sweep line x value used to compare segments'''

    OFFSET = 5
    '''*Big* offset used to compare and intersect segments'''

    EPSILON = 0.000001
    '''Threadhold to compare with 0, used in overloaded operators'''

    P = Point
    '''Alias for Point class to be used (needs to implement +, -, *, pad)'''

    def cintersect(s, o):
        "http://www.bryceboe.com/2006/10/23/line-segment-intersection-algorithm/"
        if o is None or s is None: return False

        def ccw(A, B, C): return (C.y - A.y) * (B.x - A.x) > (B.y - A.y) * (C.x - A.x)
        return ccw(s.left, o.left, o.right) != ccw(s.right, o.left, o.right) and \
           ccw(s.left, s.right, o.left) != ccw(s.left, s.right, o.right)

    def intersect(s, o):
        '''Intersect two segment intances'''
        da = s.right - s.left; da = s.P(da.y, -da.x)
        db = o.right - o.left
        return o.left + db * ((da * (s.left - o.left)) / (da * db))

    @property
    def repr(s):
        "Finds a *representant* of the segment (the cross with the sweep line)"
        My = max(abs(p.y) for p in s) + s.OFFSET
        return Segment(s.P(s.X, -My), s.P(s.X, My)).intersect(s).y

    def __str__(s): return "[%s]" % ','.join([str(s.left), str(s.right)])
    def __repr__(self): return self.__str__()

    def __cmp__(s, o):
        '''
        Compares segments s and o.
        It uses class variable `segment.X` as the point to cross the segment.

        A special attribute `tiebreak` is used for the case where 2 points
        are considered equal (intersection points). Used to decide the one
        above and below (for the tree structure to work correctly).
        '''
        if s.X is None: raise TypeError("self.X not set")
        cmp = s.repr - o.repr  # (s.repr - o.repr).y

        tiebreak = getattr(s, 'tiebreak', False)
        if tiebreak and abs(cmp) < s.EPSILON and s != o: return tiebreak
        return cmp

    def __lt__(s, o): return s.__cmp__(o) < s.EPSILON
    def __gt__(s, o): return s.__cmp__(o) > -s.EPSILON
    def __le__(s, o): return s.__cmp__(o) <= s.EPSILON
    def __ge__(s, o): return s.__cmp__(o) >= -s.EPSILON


class Event(namedtuple("Event_tuple", ['point', 'type', 'segment'])):
    ''''Event for the event queue, of type (left, right, inters)'''


class Util:
    def gkey(op, *args, **kwargs):
        "Applies args to op and filters KeyError (for tree operations)"
        try: return op(*args, **kwargs)
        except KeyError: return None

    def intersect(below, above, q, qkeys):
        "Intersects 2 segments (if possible) and inserts it in the event queue"
        if Segment.cintersect(below, above):
            intr = below.intersect(above)
            if intr not in qkeys:  # check if the point was already found
                q.put(Event(intr, "inters", [below, above]))
                qkeys.add(intr)


def segment_intersect(segments, eps=0.):
    '''
    Preconditions of segement:

    1. No two line segment endpoints or crossings have the same x-coordinate
    2. No line segment endpoint lies upon another line segment
    3. No three line segments intersect at a single point.

    Arguments
    ---------
    - `segments`: a list of segments, where each segment is represented by 2
        points (a pair of (x,y) values)

    - `eps`: the epsilon zero threashold to be used

    Return
    ------
    A list of points from intersecting the segments
    '''
    segments = [Segment(*sorted(Point(*[float(f) for f in p]) for p in s)) for s in segments]
    events = [Event(s.left, 'left', s) for s in segments] +\
             [Event(s.right, 'right', s) for s in segments]

    ilist, tree, q, qk = [], btree(), Q.PriorityQueue(), set()
    for e in events: q.put(e); qk.add(e.point)

    while not q.empty():
        point, type, seg = q.get()
        Segment.X = point.x  # set current sweep line

        if type == 'left':
            tree.insert(seg, None)
            Util.intersect(Util.gkey(tree.prev_key, seg), seg, q, qk)
            Util.intersect(seg, Util.gkey(tree.succ_key, seg), q, qk)

        elif type == 'right':
            Util.intersect(Util.gkey(tree.prev_key, seg),
                           Util.gkey(tree.succ_key, seg),
                           q, qk)
            tree.remove(seg)

        else:  # an intersection point
            ilist.append(point)
            below, above = seg  # in this case there are 2 incident segments

            # below = above, use tiebreak to reverse order
            below.tiebreak, above.tiebreak = -1, 1; tree.remove(below)

            # make bellow bigger and insert it again
            below.tiebreak, above.tiebreakm, below, above = 1, -1, above, below
            tree.insert(above, None)

            bbelow, aabove = Util.gkey(tree.prev_key, below), \
                             Util.gkey(tree.succ_key, above)
            Util.intersect(bbelow, below, q, qk)
            Util.intersect(above, aabove, q, qk)

            for a in [above, below]: delattr(a, 'tiebreak')  # undo tiebreak

    # end while not q.empty():
    return ilist


def main(args=None):
    '''
    Main function that...
    '''
    args = args or sys.argv[1:]
    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG, format="%(message)s")

    from random import shuffle

    npo = 6
    rx = list(range(2 * npo)); shuffle(rx)
    ry = list(range(2 * npo)); shuffle(ry)
    r = list(zip(rx, ry))
    rand = zip(r[0::2], r[1::2])

    rand = [Segment(*[Point(*[float(f) for f in p]) for p in s]) for s in rand]  # cast to python types

    def _plot():
        for s in rand:
            plt.plot([p[0] for p in s], [p[1] for p in s], "g")
            plt.plot(*(tuple() + s.left + ('b.', )))
            plt.plot(*(tuple() + s.right + ('r.', )))
        return rand
    global plot  # export plot for debugging purposes
    plot = _plot

    intr = segment_intersect(rand)  # ;print(intr)

    plot()
    plt.plot([p[0] for p in intr], [p[1] for p in intr], "g.")
    plt.show()

    return (rand, intr, plot)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyError:
        logger.error("KeyError, try using global plot() to debug the points in"
                     "interactive mode")
        plot()
        plt.show()
