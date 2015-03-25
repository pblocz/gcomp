#! /usr/bin/env python
# -*- encoding: utf-8 -*-

class IntersectionBezier(object):
    '''
    self.polyA
    self.polyB
    
    '''

    _EPSILON = 0.01

    def _division(self, poly, t = 0.5):
        N = len(poly)
        divA, divB = np.empty((N,2)), np.empty((N,2))
        t = np.arrary([t]); rt = 1 - t

        divA[0], divB[0] = poly[0], poly[-1]


        prevx = poly[:,0][None].T; actx = np.empty(prevx.shape)
        prevy = poly[:,1][None].T; acty = np.empty(prevy.shape)
        for j in range(1,self.N):

            np.multiply(np.roll(prevx,-1,axis=0),t,actx)
            np.add(actx,np.multiply(prevx,rt),actx)

            np.multiply(np.roll(prevy,-1,axis=0),t,acty)
            np.add(acty,np.multiply(prevy,rt),acty)

            prevx, actx = actx[:-1], prevx[:-1]
            prevy, acty = acty[:-1], prevy[:-1]

        return np.vstack((prevx[0,:],prevy[0,:])).T


    def __call__(self, polygA, polygB, epsilon = IntersectionBezier._EPSILON):
        '''
        args:
        - polygA: shape = (n,2)
        - polygB: shape = (m,2)
        
        kwargs:
        - epsilon: zero threashold

        return: k intersection points numpy array with shape = (k,2)
        '''
        pass

    def __plot__(self, k = 3, figure = None):
        '''
        Plots the curve with k steps

        kwargs:
        - k: the number of subdivision steps
        - figure: the matplotlib figure to plot
        '''
        
        pass

# Tiempos del profe
# 2.85 3.10, 3.39

import sys

def main(args=None):
    args = sys.argv
    # Generamos interactivo
    pass

if __name__ == "__main__": sys.exit(main())
