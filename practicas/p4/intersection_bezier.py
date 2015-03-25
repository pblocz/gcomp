#! /usr/bin/env python
# -*- encoding: utf-8 -*-

class IntersectionBezier(object):
    _EPSILON = 0.01

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
