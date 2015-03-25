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

    @staticmethod
    def _bbox(poly):
        return [(np.amin(poly[:,0]),np.amin(poly[:,1]),),
                (np.amax(poly[:,0]),np.amax(poly[:,1]),)]

    @staticmethod
    def _bbox_intersect(self, box1, box2):
        
        

    def __call__(self, polyA, polyB, epsilon = _EPSILON):
        '''
        args:
        - polyA: shape = (n,2)
        - polyB: shape = (m,2)
        
        kwargs:
        - epsilon: zero threashold

        return: k intersection points numpy array with shape = (k,2)
        '''
        pass

    def __plot__(self, polyA, polyB, k = 3, figure = None):
        '''
        Plots the curve with k steps

        kwargs:
        - k: the number of subdivision steps
        - figure: the matplotlib figure to plot
        '''
        if k == 0: 
            plt.plot(polyA[:,0], polyA[:,1])
            plt.plot(polyA[:,0], polyA[:,1],"ro")
            return
        
        # plt.plot(polyA[:,0], polyA[:,1])
        divA, divB = self._division(polyA)

        # print divA
        # print divB
        
        # plt.plot(divA[:,0], divA[:,1])
        # plt.plot(divB[:,0], divB[:,1])

        self.__plot__(divA, None, k - 1)
        self.__plot__(divB, None, k - 1)

        pass

# Tiempos del profe
# 2.85 3.10, 3.39

import sys

def main(args=None):
    args = sys.argv
    # Generamos interactivo

    puntos = [[0,0],[1,1],[0.2,-10],[2,0]]

    IntersectionBezier().__plot__(np.array(puntos),None, k = 4)
    plt.show()

    pass

if __name__ == "__main__": sys.exit(main())
