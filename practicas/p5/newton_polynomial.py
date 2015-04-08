#! /usr/bin/env python2
# -*- encoding: utf-8 -*-

'''
authors: Pablo Cabeza & Diego González
'''

import numpy as np
import scipy.interpolate as sc
import matplotlib.pyplot as plt

    
def divdiff(X,Y):
    '''
    Compute the divided difference of the pairs zip(X,Y)

    returns: a ndarray where r[i] = [y_0...y_i]
    '''
    res = np.empty(len(X))
    res[0] = Y[0]

    for i in range(1,len(X)):        
        Y = (Y[1:] - Y[:-1]) / (X[i:] - X[:-i])
        res[i] = Y[0]
    return res
    
    
def newton_polynomial(y, tau, num_points=100, libraries=False):
    '''
    Computes de Newton's polynomial interpolating values x at knots tau
    y: numpy array of size n; points to interpolate
    tau: numpy array of size n; knots tau[0] < tau[1] < ... < tau[n-1]
    num_points: number of points at which the polynomial will be
                evaluated

    libraries: False means only linear algebra can be used
               True means every module can be used.

    returns:
       numpy array of size num_points given by the polynomial 
       evaluated at np.linspace(tau[0], tau[1], num_points)

    Maximum cost allowed: 5,43 s at lab III computers
            degree = n - 1 = 9
            num_points = 100
    '''    
    t = np.linspace(tau[0], tau[-1], num_points)

    if libraries == False:
        dd = divdiff(tau,y)
        res = np.empty((len(tau), num_points))
        res[0].fill(1) # i = 0 -> base case

        tau_a = np.subtract(t, tau[:-1][None].T, res[1:]) # tau_a[i] = res[i+1] = t - tau[i] 
        np.multiply.accumulate(tau_a, out = tau_a) # tau_a[i] = mult[j=0,i](t - tau[j])
        np.multiply(res, dd[None].T, res) # res[i] = tau_a[i]*dd[i] = prod i factor of the newton polynomial

        return res.sum(0)

    else:
        return np.polyval(np.polyfit(tau,y,len(tau)-1), t)

if __name__ == '__main__':

    n = 10.0 # 19.0 ?!
    tau = np.arange(n)
    x = np.random.randint(-n, n, size=n)
    num_points = 100
    poly_1 = newton_polynomial(x, tau, num_points, libraries=True)
    poly_0 = newton_polynomial(x, tau, num_points, libraries=False)
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
