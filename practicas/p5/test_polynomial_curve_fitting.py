import numpy as np
import matplotlib.pyplot as plt
import timeit
from polynomial_curve_fitting import polynomial_curve_fitting


def test_least_squares_fitting():
    n = 10  
    x = np.random.randint(-10, 10, size=(n, 2))
    knots = np.linspace(0, 1, n)
    num_points = 200

    poly_0 = polynomial_curve_fitting(x, knots, method='least_squares',
                                      libraries=False, num_points=num_points)
    poly_1 = polynomial_curve_fitting(x, knots, method='newton',
                                      libraries=True, num_points=num_points)    
    t = np.linspace(knots[0], knots[-1], num_points) 

    plt.plot(poly_0[:, 0], poly_0[:, 1])
    plt.plot(poly_1[:, 0], poly_1[:, 1])
    plt.plot(x[:, 0], x[:, 1], 'o')
    plt.show()
    assert np.max(poly_0 - poly_1) < 1e-1, 'least squares fitting\
                                                 polynomial is not correct'


def test_least_squares_fitting_regularized():
    n = 10  
    x = np.array([[-7, -4], [-3, -1], [1, 2], [2, 4], [2.6, 3], [4, 1],
                  [10, 1], [12, 1], [12.4, -11], [20, -1]])
    list_L = [10**k for k in range(-15, -5)]
    print list_L
    knots = np.linspace(0, 1, n)
    num_points = 200
             
    for L in list_L:
        poly = polynomial_curve_fitting(x, knots, method='least_squares', L=L,
                                      libraries=False, num_points=num_points)
        plt.plot(poly[:, 0], poly[:, 1])
    plt.plot(x[:, 0], x[:, 1], 'o')
    
    plt.show()  
    

def test_newton_poly_cheb():
    n = 15  
    x = np.random.randint(-10, 10, size=(n, 2))
    num_points = 200
    cheb_knots = np.array([ 0.99726095,  0.97552826,  0.9330127,   0.87157241,\
                            0.79389263,  0.70336832,  0.60395585,  0.5,\
                            0.39604415,  0.29663168,  0.20610737,  0.12842759,\
                            0.0669873,   0.02447174,  0.00273905])
    poly_0 = polynomial_curve_fitting(x, 'chebyshev', method='newton',
                                      libraries=False, num_points=num_points)
    poly_1 = polynomial_curve_fitting(x, cheb_knots, method='newton',
                                      libraries=True, num_points=num_points)
    assert np.max(poly_0 - poly_1) < 1e-1, 'wrong newton polynomial with cheb_knots'
        
    plt.plot(poly_0[:, 0], poly_0[:, 1])
    plt.plot(poly_1[:, 0], poly_1[:, 1])
    plt.plot(x[:, 0], x[:, 1], 'o')
    plt.show()  
    
def test_newton_poly():
    n = 10  
    x = np.random.randint(-10, 10, size=(n, 2))
    knots = np.linspace(0, 1, n)
    num_points = 200
    
    poly_0 = polynomial_curve_fitting(x, knots, method='newton',
                                      libraries=False, num_points=num_points)
    poly_1 = polynomial_curve_fitting(x, knots, method='newton',
                                      libraries=True, num_points=num_points)
    assert np.linalg.norm(poly_0 - poly_1) < 1e-2, 'wrong newton polynomial'
        
    
    plt.figure()
    plt.plot(poly_0[:, 0], poly_0[:, 1])
    plt.plot(poly_1[:, 0], poly_1[:, 1])
    plt.plot(x[:, 0], x[:, 1], 'o')
    
    plt.show()  

def timing_curve_fitting():
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

if __name__ == '__main__':
    
    test_least_squares_fitting()
    test_least_squares_fitting_regularized()
    #test_newton_poly_cheb()
    test_newton_poly()
    timing_curve_fitting()
