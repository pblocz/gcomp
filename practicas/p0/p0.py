# -*- coding: utf-8 -*-
'''
Practica 0

@author: Pablo Cabeza García
'''

from scipy.integrate import odeint # instalarla
import matplotlib.pyplot as plt
import numpy as np
from math import *

def solve(k,I,s0,delta=0.1,init_cond=(0,0,1,0)):
    '''
    @arg k: funcion de curvatura
    @arg I: intervalo [x0,xN]
    @arg s0:  punto donde se dan las cond iniciales

    @kwarg delta: para el intervalo
    @kwarg init_cond: [x,y,dx,dy] iniciales
    '''

    def rhs_eqs(Y,s):
        x, y , dx, dy = Y
        return [dx, dy, -k(s) * dy, k(s) * dx] 

    I0,In = I
    intervs = [np.arange(s0,I0-delta,-delta), np.arange(s0,In + delta, delta)]
    solvs = [odeint(rhs_eqs, init_cond, interv) for interv in intervs]
    solv = np.concatenate((solvs[0][::-1],solvs[1]),axis=0)
    

    print intervs
    
    curve_x = solv[:,0]; curve_y = solv[:,1]

    plt.axes().set_aspect('equal')
    plt.plot(curve_x, curve_y)
    plt.show()

if __name__ == "__main__":
    # solve(lambda s: 1, [0, 2*np.pi], np.pi)
    # solve(lambda s: s * (s - 3) * (s + 5), [0,10], 0, delta=0.01)
    solve(lambda s: exp(s), [-2*np.pi,2*np.pi], 0)

    #solve(lambda s: 4*cos(s)**2 + sin(s)**2, [-100,100], 0)
