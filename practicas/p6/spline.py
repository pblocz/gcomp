def spline2d(a, b, psi, k, nu, A, num_dots):
    '''
    Computes a plane spline curve of order k defined on the interval 
    [a, b] with knots psi, multiplicities nu and coefficients A.

    Parameters:
      - a, b -- ends of the interval, real numbers
      - psi -- list of breakpoints, a < psi[0] < ... < psi[-1] < b
      - k -- order of the curve, the degree is <= k - 1
      - nu -- list of integer multiplicities of each breakpoint,
    	len(psi) = len(nu), 1 <= nu[i] < k
      - A -- list of coefficients of the B-spline basis,
    	A = [[x0, y0], [x1, y1],..., [xN, yN]
      - num_dots -- number of dots of the spline to be plotted,
    	uniformly spaced along the interval [a, b]

    Returns:
    the spline curve as a numpy array of size (2, num_dots) 
    '''

    t = np.linspace(a,b,num_dots)
