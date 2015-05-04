import numpy as np
import matplotlib.pyplot as plt

IF = 0

def boor_coef(i, r, k, knots, t):
    global IF
    # print IF*'\t', "w:: i, r, k::", i,r, k
    denom = float(knots[i + k + 1 - r] - knots[i])
    if denom == 0: return 0
    return (t - knots[i]) / denom


def boor_step(i, r, k, knots, points, t):
    global IF
    # print IF*'\t', "i, r, k::", i, r, k
    if r == 0:
        ret = np.empty_like(t); ret.fill(points[i])
        # print IF*'\t',"result::", ret
        return ret

    IF+=1
    coef = boor_coef(i, r, k, knots, t)
    IF-=1
    # print coef
    IF += 1
    boor1 = boor_step(i - 1, r - 1, k, knots, points, t)
    boor2 = boor_step(i, r - 1, k, knots, points, t)
    IF -= 1
    res = (1 - coef) * boor1 +\
            coef * boor2
    # print IF*'\t',"result::", res
    return res



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

    t = np.linspace(a, b, num_dots)

    nknots = 2*k + k*len(nu) - sum(nu)
    knots, k_points = np.empty(nknots), np.empty([nknots, 2])

    knots[:k].fill(a); knots[-k:].fill(b)
    k_points[:k][:, 0].fill(A[0, 0]); k_points[:k][:, 1].fill(A[0, 1])
    k_points[-k:][:, 0].fill(A[-1, 0]); k_points[-k:][:, 1].fill(A[-1, 1])

    idx = k
    for i, n in enumerate(nu):
        off = k - n
        knots[idx:idx + off].fill(psi[i])
        k_points[idx:idx + off][:, 0].fill(A[i, 0])
        k_points[idx:idx + off][:, 1].fill(A[i, 1])
        idx += off


    # print "m, n, k", len(knots), len(A), k, "::", len(knots), "==", len(A)+k+1

    ret1, ret2, N = [], [], len(knots)
    for l in xrange(N):  # calculamos a_j^{k-1}
        if l == N - 1: break


        cond = (t >= knots[l]) & (t < knots[l+1]) #\
            # ((t <= knots[l+1]) if l == N - 2 else (t < knots[l+1]))
        eval_points = t[cond]

        # print "eval_points::", eval_points
        # print "l, knots[l], knots[l+1]::", l, knots[l], knots[l+1]
        if len(eval_points) == 0: continue

        IF=0
        ret1.append(boor_step(l, k - 1, k - 1, knots, A[:,0], eval_points))
        IF=0
        ret2.append(boor_step(l, k - 1, k - 1, knots, A[:,1], eval_points))

    print t
    print [np.concatenate(ret1), np.concatenate(ret2)]
    return np.array([np.concatenate(ret1), np.concatenate(ret2)]).T

if __name__ == "__main__":
    '''
    num_points = 100
    k = 4
    A = np.array([[0,0], [1,1],[2,2],[3,2],[4,2],[5,1],[6,0],[7,1]]).astype(float)
    a, b, = 0., 1.
    xi = np.array([0.25, 0.5, 0.75])
    nu = np.array([2, 2, 2]).astype(float)  # [3,3,3])
    '''

    a = 0
    b = 4
    xi = np.array([1, 2, 3])
    nu = np.array([2, 2, 2])
    k = 3
    A = np.array([[-3, 3], [-3, 3], [-3, 3], [8, 2], [-4, 8], [0, 8], [-1, 0], [-1, 0], [-1, 0]])
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

    r = spline2d(a, b, xi, k, nu, A, num_points)
    plt.plot(A[:,0], A[:,1], 'ro-')
    plt.plot(r[:,0],r[:,1], 'b')
    plt.show()
