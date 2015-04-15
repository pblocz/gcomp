from __future__ import division
import numpy as np
import matplotlib.pyplot as plt



def polynomial_curve_fitting(points, knots, method, L=0, libraries=False,
                             num_points=100):    
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
    pass


def polynomial_curve_fitting1d(points, knots, method, L=0, libraries=False,
                             num_points=100):    
    pass


def newton_polynomial(x, tau, num_points=100, libraries=False):    
    #I've used polyfit and polyval if libraries==True
    #If you find something faster, I'll be glad to know about it.
    pass

def eval_poly(t, coefs, tau=None):    
    pass
        
def least_squares_fitting(points, knots, num_points, L=0, libraries=True):    
    #I've used np.linalg.lstsq and np.polyval if libraries==True
    pass
        
def chebyshev_knots(a, b, n):
    pass
       



                                 
                
