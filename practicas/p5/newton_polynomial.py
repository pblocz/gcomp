#! /usr/bin/env python2
# -*- encoding: utf-8 -*-

'''
authors: Pablo Cabeza & Diego González
'''

import numpy as np
import scipy.interpolate as sc
import matplotlib.pyplot as plt


def polynomial_curve_fitting(points, knots, method, *args, **kwargs): # L=0, libraries=False, num_points=100):    
    '''
       Fits planar curve to points at given knots. 

       Arguments:
           points -- coordinates of points to adjust (x_i, y_i) given by a numpy array of shape (N, 2)
           knots -- strictly increasing sequence at which the curve will fit the points, tau_i
               It is given by a np.array of shape M, unless knots='chebyshev', in this case
                   N Chebyshev's nodes between 0 and 1 will be used instead of tau.
           method -- one of the following: 
               'newton' computes the interpolating polynomial curve using Newton's method.
                   returns error if N!=M. 
               'least_squares' computes the best adjusting curve in the least square sense,
                   i.e., min_a ||Ca - b||**2 + L/2 ||a||**2
           L -- regularization parameter
           libraries -- If False, only numpy linear algebra operations are allowed. 
               If True, any module can be used. In this case, a very short and fast code is expected
           num_points -- number of points to plot between tau[0] and tau[-1]

       Returns:
           numpy array of shape (num_points, 2) given by the evaluation of the polynomial
           at the evenly spaced num_points between tau[0] and tau[-1]
    '''
    if method == "newton": return newton_polynomial(points, knots, *args, **kwargs)
    
def divdiff(X,Y):
    '''
    Compute the divided difference of the pairs zip(X,Y)

    returns: a ndarray where r[i] = [y_0...y_i]
    '''
    res = np.empty(len(X))
    res[0] = Y[0]


    for i in xrange(1,len(X)):        
        Y = (Y[1:] - Y[:-1]) / (X[i:] - X[:-i])
        res[i] = Y[0]
    return res

def newton_polynomial1d(y,tau, num_points, out = None):
    '''
    It uses naïve polynomial evalutation instead of horner algorithm. For
    each t it does 2*n-1 multiplications, n-1 additions and n-1
    subtracttions, which translates into 5 numpy calls with no python 
    for-loops.

    It runs in (except for divdiff()):

        O(n*p) operations inside numpy
        O(1) numpy calls
        O(n*p) space storage 
    '''

    t = np.linspace(tau[0], tau[-1], num_points)
    if out is None: out = np.empty(num_points)

    dd = divdiff(tau,y)
    res = np.empty((len(tau), num_points))

    res[0].fill(1) # i = 0 -> base case

    tau_a = np.subtract(t, tau[:-1][None].T, res[1:]) # tau_a[i] = res[i+1] = t - tau[i] 
    np.multiply.accumulate(tau_a, out = tau_a) # tau_a[i] = mult[j=0,i](t - tau[j])
    np.multiply(res, dd[None].T, res) # res[i] = tau_a[i]*dd[i] = factor i of the newton polynomial

    return res.sum(0, out = out) # sum of factors

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
        '''
        It uses naïve polynomial evalutation instead of horner algorithm. For
        each t it does 2*n-1 multiplications, n-1 additions and n-1
        subtracttions, which translates into 5 numpy calls with no python 
        for-loops.

        It runs in (except for divdiff()):

            O(n*p) operations inside numpy
            O(1) numpy calls
            O(n*p) space storage 
        '''

        r_val = np.empty((2, num_points))
        newton_polynomial1d(y[:,0],tau, num_points, out = r_val[0])
        newton_polynomial1d(y[:,1],tau, num_points, out = r_val[1])
        return r_val.T
        
    else: 
        fit = np.polyfit(tau,y,len(tau)-1)
        return np.vstack([np.polyval(fit[:,i],t) for i in [0,1]]).T


if __name__ == '__main__':

    n = 10.0 # 19.0 ?!
    tau = np.arange(n)
    x = np.random.randint(-n, n, size=(n,2))
    num_points = 100
    poly_1 = newton_polynomial(x, tau, num_points, libraries=True)
    poly_0 = newton_polynomial(x, tau, num_points, libraries=False)
    print np.linalg.norm(poly_0 - poly_1)
    
    t = np.linspace(tau[0], tau[-1], num_points)    
    plt.plot(poly_0[:,0],poly_0[:,1], 'b')
    plt.plot(poly_1[:,0],poly_1[:,1], 'r')
    plt.plot(x[:,0],x[:,1], 'o')

    plt.show()
    
    
    import timeit
    
    print(timeit.repeat("x = np.random.randint(-10, 10, size=n); newton_polynomial(x, tau, libraries=False)",
                        setup="from __main__ import newton_polynomial, n,  tau, np",
                        number=10000))
    print(timeit.repeat("x = np.random.randint(-10, 10, size=n); newton_polynomial(x, tau, libraries=True)",
                        setup="from __main__ import newton_polynomial, n,  tau, np",
                        number=10000))
