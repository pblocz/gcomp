#! env python2
# -*- encoding: utf-8 -*-

'''
authors: Pablo Cabeza & Diego González
'''

import numpy as np
import scipy.interpolate as sc
import matplotlib.pyplot as plt

def dif_div(X,Y):
    res = np.empty(len(X))
    res[0] = Y[0]
    X_u = X
    for i in range(len(X)-1):
        X_u = np.roll(X_u,-1)
        Y = (np.roll(Y,-1) - Y) / (X_u - X)
        res[i+1] = Y[0]
    return res
    

def newton_polynomial(y, tau, num_points=50, libraries=False):
    dd = dif_div(tau,y)
    t = np.linspace(tau[0], tau[-1], num_points)    
    
    #your code here
    if libraries == False:
        res = np.empty(num_points); res.fill(dd[0])
        tau_a = np.ones(num_points)
        for i in range(1,len(tau)):
            tau_a = tau_a * (t - tau[i-1])
            res += tau_a * dd[i]
            
        return res #np.array of size num_points

    else:
        interp = sc.interp1d(tau,y, kind=len(tau)-1)
        return interp(t) #np.array of size num_points
    

import sys
if __name__ == '__main__':
    '''
    p0 = np.array([0.,1.,2.,3.])
    py = np.array([2.,4.,6.,1.])

    poly_0 = newton_polynomial(py,p0,50, libraries=True)
    print poly_0
    t = np.linspace(p0[0], p0[-1], 50)    
    plt.plot(t, poly_0)
    plt.plot(p0, py, 'o')
    plt.show()
    sys.exit(0)
    '''
    n = 10.0
    tau = np.arange(n)
    x = np.random.randint(-n, n, size=n)
    num_points = 100
    poly_0 = newton_polynomial(x, tau, num_points, libraries=False)
    poly_1 = newton_polynomial(x, tau, num_points, libraries=True)
    print np.linalg.norm(poly_0 - poly_1)
    
    t = np.linspace(tau[0], tau[-1], num_points)    
    plt.plot(t, poly_0, 'b')
    plt.plot(t, poly_1, 'r')
    plt.plot(tau, x, 'o')

    plt.show()
    
    
    import timeit
    
    print(timeit.repeat("x = np.random.randint(-10, 10, size=n); newton_polynomial(x, tau, libraries=False)",
                        setup="from __main__ import newton_polynomial, n,  tau, np",
                        number=10000))
    print(timeit.repeat("x = np.random.randint(-10, 10, size=n); newton_polynomial(x, tau, libraries=True)",
                        setup="from __main__ import newton_polynomial, n,  tau, np",
                        number=10000))


# Mario@jinapp.com
