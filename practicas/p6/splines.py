#! /usr/bin/env python2
# -*- encoding: utf-8 -*-

'''
Different implementations for interpolating polynomials.

authors: Pablo Cabeza & Diego González
date: 29/04/2015 
'''

import numpy as np
import scipy.interpolate as sc
import matplotlib.pyplot as plt
import numpy.linalg as lalg


def w(i,k,t,xi):
    print "w", i,k
    if (xi[i]==xi[i+k-1]):
        return 0
    else:
        return (t-xi[i]) / (xi[i+k-1]-xi[i])

def ra(i,r,knots,k,A,t):
    "Calcula a_i^{r}"
    r = r-1; wik = w(i,k-r,t,knots)
    # print wik, i,k-r
    
    if r==0: 
        return A[i]
    else: 
        ra1 = ra(i-1,r,knots,k,A,t)
        ra2 = ra(i,r,knots,k,A,t)
        ret = (1-wik)*ra1 \
            + wik*ra2
        print "ra1,ra2", ra1,ra2
        print wik, ret
        return ret

def spline2d(a, b, xi, k, nu, A, num_dots):
    '''
    Computes a plane spline curve of order k defined on the interval [a, b] with knots psi,
    multiplicities nu and coefficients A.
    Parameters:
    a, b     --  ends of the interval, real numbers.
    xi       --  list of breakpoints, a < xi[0] < ... < xi[-1] < b.
    k        --  order of the curve, the degree is <= k - 1.
    nu       --  number of smoothness conditions satisfied by the curve at each breakpoint.
    A        --  list of coefficients of the B-spline basis, A = [[x0, y0], [x1, y1],..., [x[N], y[N]]
    num_dots --  number of dots of the spline to be plotted, uniformly spaced along the interval [a, b].
    Returns:     the spline curve as a numpy array of size (2, num_dots) 
    '''
    
    t = np.linspace(a,b,num_dots)
    '''
    nknots = 2*k + k*len(nu) - sum(nu)
    eknots = []
    knots = np.empty(nknots)
    aknots = np.empty([nknots,2])
    
    knots[:k].fill(a)
    aknots[:k][:,0].fill(A[0,0])
    aknots[:k][:,1].fill(A[0,1])
    knots[-k:].fill(b)
    aknots[-k:][:,0].fill(A[-1,0])
    aknots[-k:][:,1].fill(A[-1,1])
    
    idx = k
    for i,n in enumerate(nu):
         eknots.append(idx-1)
         knots[idx:idx + k - n].fill(xi[i])
         aknots[idx:idx + k - n][:,0].fill(A[i,0])
         aknots[idx:idx + k - n][:,1].fill(A[i,1])
         idx += k - n
    eknots.append(idx-1)

    '''
    ret1 = []
    ret2 = []
    
    knots = np.empty(len(xi)+2); 
    knots[0], knots[-1]=a,b
    knots[1:-1] = xi
    
    for j in xrange(len(knots)):
        # calculamos a_j^{k-1}
        puntos_evaluar = t[(t >= knots[j]) & (t < knots[j+1])] 
        ret1.append(ra(j,k-1,knots,k,A[:,0],puntos_evaluar))
        ret2.append(ra(j,k-1,knots,k,A[:,1],puntos_evaluar))
    
    return np.array([np.concatenate(ret1), np.concatenate(ret2)])
        
   
    


if __name__ == '__main__':

    a, b = 0,1
    k = 3
    xi = np.linspace(a+0.1,b-0.1,k)
    nu = np.random.randint(1,k,size=k)
    A = np.random.randint(-10, 10, size=(k+2,2))
    num_points = 100

    '''
    poly_1 = polynomial_curve_fitting(x, tau, 'least_squares', num_points, degree = None, L=0, libraries = False)
    poly_0 = polynomial_curve_fitting(x, tau, 'newton', num_points, libraries=False, L=0, degree = 2)
    print np.linalg.norm(poly_0 - poly_1)
    
    t = np.linspace(tau[0], tau[-1], num_points)  
    plt.plot(poly_0[:,0],poly_0[:,1], 'b')
    plt.plot(poly_1[:,0],poly_1[:,1], 'r')
    plt.plot(x[:,0],x[:,1], 'o')
    '''
    
    
    r = spline2d(a,b,xi,k, nu,A,num_points)
    plt.plot(A[:,0], A[:,1], 'ro-')
    plt.plot(r[:,0],r[:,1], 'b')
    plt.show()