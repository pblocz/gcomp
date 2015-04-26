#! /usr/bin/env python2
# -*- encoding: utf-8 -*-

'''
Different implementations for interpolating polynomials.

authors: Pablo Cabeza & Diego González
date: 20/04/2015 
'''

import numpy as np
import scipy.interpolate as sc
import matplotlib.pyplot as plt
import numpy.linalg as lalg

def polynomial_curve_fitting(points, knots, method, num_points=100, libraries=False, L=0, degree=None): # L=0, libraries=False, num_points=100):    
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
           degree -- degree of the polynomial. Needed only if method='least_squares'.
               If degree=None, the function will return the interpolating polynomial.


       Returns:
           numpy array of shape (num_points, 2) given by the evaluation of the polynomial
           at the evenly spaced num_points between tau[0] and tau[-1]
    '''

    if knots == 'chebyshev': knots = chebyshev_knots(0,1,len(points))

    if method == "newton": return newton_polynomial(points, knots, num_points, libraries)
    elif method == "least_squares": 
      return least_squares_fitting(points, knots, degree or (len(knots) - 1), num_points,  L, libraries)


def divdifnd(y,tau):
  "Compute divided differences for y of dimension (n,m), m > 1"

  dd = np.empty(y.shape); dd[0] = y[0]
  X = np.tile(tau,(y.shape[1],1,)).T
  yaux = y

  for i in xrange(1,len(tau)):        
      yaux = (yaux[1:] - yaux[:-1]) / (X[i:] - X[:-i])
      dd[i] = yaux[0]

  return dd


def newton_polynomialnd(y,tau, num_points, out = None):
    '''
    It uses naïve polynomial evalutation instead of horner algorithm. For
    each t it does 2*n-1 multiplications, n-1 additions and n-1
    subtracttions, which translates into 5 numpy calls with no python 
    for-loops.

    arguments:
      - out: numpy.array of dimension (num_points, y.shape[1]) to fill or none

    return: a reference to out if given, or a new numpy.array with the result

    It runs in (except for divdifnd()):

        O(n*p) operations inside numpy
        O(1) numpy calls
        O(n*p) space storage 
    '''

    t = np.linspace(tau[0], tau[-1], num_points)
    if out is None: out = np.empty((num_points,y.shape[1],))

    dd = divdifnd(y,tau)

    # prepare `res` to hold `res[i+1] = (t - tau[0])*...*(t-tau[i])`
    res = np.empty((len(tau), num_points)); res[0].fill(1) # i = 0 -> base case
    tau_a = np.subtract(t, tau[:-1][None].T, res[1:]) # tau_a[i] = res[i+1] = t - tau[i] 
    np.multiply.accumulate(tau_a, out = tau_a) # tau_a[i] = mult[j=0,i](t - tau[j])
    
    # For each column of points, use precomputed res to get the evaluation in `t`
    aux = np.empty_like(res)
    for i in xrange(dd.shape[1]):
      np.multiply(res, dd[:,i][None].T, out = aux) # res[i] = tau_a[i]*dd[i] = factor i of the newton polynomial
      aux.sum(0, out = out[:,i]) # sum of factors

    return out

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
    '''

    if libraries == False: return newton_polynomialnd(y, tau, num_points)    
    else: 
      t = np.linspace(tau[0], tau[-1], num_points)
      fit = np.polyfit(tau,y,len(tau)-1)
      return np.vstack([np.polyval(fit[:,i],t) for i in xrange(fit.shape[1])]).T


def least_squares_fitting(points, knots, degree, num_points, L=0, libraries=True):
    '''
    Computes the interpolating curve given by the polynomial of the given degree.

    arguments:
      - points: numpy array of size n; points to interpolate
      - knots: ordered numpy array of size n
      - num_points: number of points at which the polynomial will be evaluated
      - degree: the degree of the polynomial

    keyword arguments:
      - libraries: False means only linear algebra can be used
               True means every module can be used.
      - L: the correction parameter to avoid overfeeting in the polynomial

    returns:
       numpy array of size num_points given by the polynomial 
       evaluated at np.linspace(knots[0], knots[-1], num_points)
    '''

    t = np.linspace(knots[0], knots[-1], num_points)

    if libraries == False:
      # Compute vardemond matrix transposed `C` (fast) (in the algo this is `C'`)
      C = np.tile(knots,(degree+1,1)); C[0].fill(1)
      np.multiply.accumulate(C[1:], out = C[1:])

      # Compute `H` matrix and polynomial coefficients `a`
      H = np.dot(C,C.T)+ L/2.0 * np.identity(degree+1)
      a = np.dot(np.dot(lalg.inv(H),C),points)

      # Eval polynomial in `t`, using fast numpy operations and O(n*p) memory
      T = np.tile(t,(degree+1,1)); T[0].fill(1)
      np.multiply.accumulate(T[1:], out = T[1:])
      return np.dot(T.T, a)

    else: 
      C = np.vander(knots, int(degree) + 1)[:,::-1] # Slower than our computing, left because is clearer

      H = np.dot(C.T,C)+ L/2.0 * np.identity(degree+1)
      fit = np.linalg.lstsq(H,np.dot(C.T,points))[0][::-1]

      return np.vstack([np.polyval(fit[:,i],t) for i in [0,1]]).T


def chebyshev_knots(a, b, n):
  '''
  Computes `n` chebyshev_knots in the interval `[a,b]`

  returns: numpy.array of the nodes
  '''
  return ((a+b) + (b-a) * np.cos( (2*np.arange(1,n+1) - 1)*np.pi / (2.0*n) )) / 2.0



if __name__ == '__main__':

    n = 10.0
    tau = np.arange(1,n+1)
    x = np.vstack((tau,np.random.randint(-n, n, size=n))).T
    num_points = 100

    poly_1 = polynomial_curve_fitting(x, tau, 'least_squares', num_points, degree = None, L=0, libraries = False)
    poly_0 = polynomial_curve_fitting(x, tau, 'newton', num_points, libraries=False, L=0, degree = 2)
    print np.linalg.norm(poly_0 - poly_1)
    
    t = np.linspace(tau[0], tau[-1], num_points)  
    plt.plot(poly_0[:,0],poly_0[:,1], 'b')
    plt.plot(poly_1[:,0],poly_1[:,1], 'r')
    plt.plot(x[:,0],x[:,1], 'o')

    plt.show()
    
    import timeit

    setup = '''
import numpy as np
from polynomial_curve_fitting import polynomial_curve_fitting
n=10
knots = np.linspace(0, 1, n)
num_points = 200
    '''
    print 'newton:', min(timeit.repeat("x = np.random.randint(-10, 10, size=(n, 2));\
        polynomial_curve_fitting(x, knots, method='newton',\
        libraries=False, num_points=num_points)", setup=setup, number=10000))
        
    print 'newton_libraries:', min(timeit.repeat("x = np.random.randint(-10, 10, size=(n, 2));\
        polynomial_curve_fitting(x, knots, method='newton',\
        libraries=False, num_points=num_points)", setup=setup, number=10000))
        
    print 'least_squares:', min(timeit.repeat("x = np.random.randint(-10, 10, size=(n, 2));\
        polynomial_curve_fitting(x, knots, method='least_squares',\
        libraries=False, num_points=num_points, L=1e-10)", setup=setup, number=10000))

    print 'least_squares_libraries:', min(timeit.repeat("x = np.random.randint(-10, 10, size=(n, 2));\
        polynomial_curve_fitting(x, knots, method='least_squares',\
        libraries=True, num_points=num_points, L=1e-10)", setup=setup, number=10000))