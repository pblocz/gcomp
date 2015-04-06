#! /usr/bin/env python
# -*- encoding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple, deque




class FastIntersectionBezierCore(object):
    '''
    Utility class that holds implementations details of the algorithm
    '''

    _EPSILON = 0.01
    
    @staticmethod
    def intersect_bbox_check(pA, pB, epsilon = _EPSILON):
        AX, AY = np.amax(pA,axis=0); Ax, Ay = np.amin(pA, axis=0)
        BX, BY = np.amax(pB,axis=0); Bx, By =  np.amin(pB, axis=0)
        return Ax <= BX and AX >= Bx and Ay <= BY and AY >= By

    @staticmethod
    def division(poly, t = 0.5):
        divA, divB = np.empty(poly.shape), np.empty(poly.shape)
        divA[0], divB[0] = poly[0], poly[-1]

        curr = np.copy(poly)
        for j in xrange(1,len(poly)): 
            np.add(curr,np.roll(curr,-1,axis=0), curr)
            curr = np.multiply(curr,t)[:-1]
            # np.multiply(curr,t,curr) # curr *= t

            divA[j], divB[j] = curr[0], curr[-1]

        return (divA, divB)


    @staticmethod
    def seg_intersect(a1,a2, b1,b2) :
        da = a2-a1; da[0], da[1] = da[1], -da[0]
        db = b2-b1
        return (np.dot(da,a1-b1)/np.dot(da,db))*db + b1
        # return db*np.dot(da,a1-b1)/np.dot(da,db) + b1


    @staticmethod
    def _delta2(poly):
        a =  np.roll(poly,-1,0)
        np.subtract(poly,a,a)

        return np.subtract(a,np.roll(a,-1,0),a)[:-2]




class IntersectionBezier(FastIntersectionBezierCore):
    '''
    self.polyA
    self.polyB
    
    '''

    _EPSILON = FastIntersectionBezierCore._EPSILON
    
    
    def __call__(self, polyA, polyB, epsilon = _EPSILON, firstCall = True):
        '''
        Computes the intersection of the two bezier curves defined by the
        given polygons. It stores the results in:
           - self.intersectionPoints
           - self.polyA
           - self.polyB

        args:
        - polyA: shape = (n,2)
        - polyB: shape = (m,2)
        
        kwargs:
        - epsilon: zero threashold

        return: k intersection points numpy array with shape = (k,2)
        '''

        if firstCall: self.lastPolyA, self.lastPolyB = polyA, polyB

        ret = []
        if self.intersect_bbox_check(polyA,polyB):
            deltaA = self._delta2(polyA)
            deltaB = self._delta2(polyB)

            DA, DB = np.amax(np.absolute(deltaA)), np.amax(np.absolute(deltaB))
            n, m, eps = len(polyA), len(polyB), epsilon*8

            if n*(n-1)*DA > eps: polyA, (divA,divB) = polyB, self.division(polyA)
            elif m*(m-1)*DB > eps: polyA, (divA,divB) = polyA, self.division(polyB)
            else: ret = [self.seg_intersect(polyA[0], polyA[-1],polyB[0], polyB[-1])] 

            if not ret: 
                ret = self(polyA,divA, epsilon = epsilon, firstCall = False)
                ret.extend(self(polyA,divB, epsilon = epsilon, firstCall = False))


        self.lastIntersectionPoints = ret
        return ret



    def __plot_points(self, polyA, k = 3):
        '''
        Compute the points of the bezier curve with k steps

        kwargs:
        - k: the number of subdivision steps
        - figure: the matplotlib figure to plot
        '''
        if len(polyA) == 0: return []
        if k == 0: 
            return polyA
        
        divA, divB = self.division(polyA)

        r = self._plot(divA, k - 1)
        return np.concatenate((r,self._plot(divB[::-1], k - 1)))

    def _plot(self, polyA, k = 3, **kwargs):
        '''
        Plots the curve with k steps

        kwargs:
        - k: the number of subdivision steps
        - figure: the matplotlib figure to plot
        '''
        if k == 0: 
            plt.plot(polyA[:,0], polyA[:,1], **kwargs)
            return
        
        divA, divB = self.division(polyA)

        self._plot(divA,  k - 1, **kwargs)
        self._plot(divB,  k - 1, **kwargs)


    def plot(self, k = 3):
        ''' Plots the last intersection curves'''
        self._plot(self.lastPolyB, k, color="b")
        self._plot(self.lastPolyA, k, color="g")

        intr = self.lastIntersectionPoints
        plt.plot([x for x,_ in intr], [y for _,y in intr], "ro")




import sys

def main(args=None):
    args = args or sys.argv
    # Generamos interactivo

    puntosA = [[0,0],[1,2],[3,1],[0,1]]
    puntosB = [[0,2],[1,0],[2,0],[0,2]]

    
    #IntersectionBezier()._plot(np.array(puntosA), k = 4)
    #IntersectionBezier()._plot(np.array(puntosB), k = 4)
    
    intersection = IntersectionBezier()
    intr = intersection(np.array(puntosA),np.array(puntosB))
    # intr = IntersectionBezier().multi_call(np.array(puntosA),np.array(puntosB))

    intersection.plot(k = 4)

    # plt.plot([x for x,_ in intr], [y for _,y in intr], "ro")
    plt.show()

if __name__ == "__main__": sys.exit(main())
