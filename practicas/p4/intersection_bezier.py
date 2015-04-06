#! /usr/bin/env python
# -*- encoding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

class IntersectionBezier(object):
    '''
    self.polyA
    self.polyB
    
    '''

    _EPSILON = 0.01
    
    from collections import namedtuple
    Rectangle = namedtuple('Rectangle', 'xmin ymin xmax ymax')
    
    def _division(self, poly, t = 0.5):
        N = len(poly)
        divA, divB = np.empty((N,2)), np.empty((N,2))
        t = np.array([t]); rt = 1 - t

        divA[0], divB[0] = poly[0], poly[-1]


        prevx = poly[:,0][None].T; actx = np.empty(prevx.shape)
        prevy = poly[:,1][None].T; acty = np.empty(prevy.shape)

        for j in range(1,N):

            actx = np.multiply(np.roll(prevx,-1,axis=0),t)
            np.add(actx,np.multiply(prevx,rt),actx)

            acty = np.multiply(np.roll(prevy,-1,axis=0),t)
            np.add(acty,np.multiply(prevy,rt),acty)

            prevx, actx = actx[:-1], prevx[:-1]
            prevy, acty = acty[:-1], prevy[:-1]

            divA[j], divB[j] = np.array([prevx[0][0],prevy[0][0]]), np.array([prevx[-1][0],prevy[-1][0]])

        return (divA, divB[::-1])

    @classmethod
    def _bbox(cls,poly):
        return cls.Rectangle(np.amin(poly[:,0]),np.amin(poly[:,1]),np.amax(poly[:,0]),np.amax(poly[:,1]))

    @staticmethod
    def _bbox_intersect(box1, box2):
        dx = min(box1.xmax, box2.xmax) - max(box1.xmin, box2.xmin)
        dy = min(box1.ymax, box2.ymax) - max(box1.ymin, box2.ymin)
        return (dx >= 0) and (dy >= 0)

    @staticmethod
    def _seg_intersect(a1,a2, b1,b2) :
        def perp( a ) :
            b = np.empty_like(a)
            b[0] = -a[1]
            b[1] = a[0]
            return b
        da = a2-a1
        db = b2-b1
        dp = a1-b1
        dap = perp(da)
        denom = np.dot( dap, db)
        num = np.dot( dap, dp )
        return (num/denom)*db + b1
    
    def __call__(self, polyA, polyB, epsilon = _EPSILON):
        '''
        args:
        - polyA: shape = (n,2)
        - polyB: shape = (m,2)
        
        kwargs:
        - epsilon: zero threashold

        return: k intersection points numpy array with shape = (k,2)
        '''
        ret = []
        
        bboxA, bboxB = self._bbox(polyA), self._bbox(polyB)
        def delta(poly):
            return (poly - np.roll(poly,-1,axis=0))[:-1]
        
        deltaA, deltaB = delta(delta(polyA)), delta(delta(polyB))
        #DA = np.apply_along_axis(np.linalg.norm, 1, deltaA)
        DA, DB = np.amax(np.sum(np.abs(deltaA)**2,axis=-1)**(1./2)), np.amax(np.sum(np.abs(deltaB)**2,axis=-1)**(1./2))
        n, m = len(polyA), len(polyB)
        if self._bbox_intersect(bboxA, bboxB):
            if n*(n-1)/160.*DA > epsilon:
                divA, divB = self._division(polyA)
                ret.extend(self(divA,polyB, epsilon = epsilon))
                ret.extend(self(divB,polyB, epsilon = epsilon))
            elif m*(m-1)/160.*DB > epsilon:
                divA, divB = self._division(polyB)
                ret.extend(self(polyA,divA, epsilon = epsilon))
                ret.extend(self(polyA,divB, epsilon = epsilon))
            else:
                return [self._seg_intersect(polyA[0], polyA[-1],polyB[0], polyB[-1])]
        return ret

    def __plot__(self, polyA, polyB, k = 3, figure = None):
        '''
        Plots the curve with k steps

        kwargs:
        - k: the number of subdivision steps
        - figure: the matplotlib figure to plot
        '''
        if k == 0: 
            plt.plot(polyA[:,0], polyA[:,1])
            #plt.plot(polyA[:,0], polyA[:,1],"ro")
            return
        
        # plt.plot(polyA[:,0], polyA[:,1])
        divA, divB = self._division(polyA)

        # print divA
        # print divB
        
        # plt.plot(divA[:,0], divA[:,1])
        # plt.plot(divB[:,0], divB[:,1])

        self.__plot__(divA, None, k - 1)
        self.__plot__(divB, None, k - 1)

        #pass

# Tiempos del profe
# 2.85 3.10, 3.39

import sys

def main(args=None):
    args = args or sys.argv
    # Generamos interactivo

    puntosA = [[0,0],[1,2],[3,1],[0,1]]
    puntosB = [[0,2],[1,0],[2,0],[0,2]]

    IntersectionBezier().__plot__(np.array(puntosA),None, k = 4)
    IntersectionBezier().__plot__(np.array(puntosB),None, k = 4)
    
    print type(np.array(puntosB))
    intr = IntersectionBezier()(np.array(puntosA),np.array(puntosB))

    plt.plot([x for x,_ in intr], [y for _,y in intr], "ro")
    plt.show()
    pass

if __name__ == "__main__": sys.exit(main())
