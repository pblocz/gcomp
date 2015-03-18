import numpy as np
from scipy.special import binom

class BezierPolicy(object):
    '''
    Base class used to implement different curve
    methods
    '''

    def __init__(self, polygon):
        self.polygon = np.array(polygon)
        self.N = len(self.polygon)

    def __call__(self, npoints=100, *args, **kwargs): 
        return self.compute(np.linspace(0, 1, npoints), 
                            *args, **kwargs)
        

class BernsteinPolicy(BezierPolicy):

    def _berstein(self,t):
        res = np.zeros((self.N+1, len(t),))
        for i in range(self.N):        
            res[i, :] = binom(self.N - 1 , i) * t**i * (1 - t)**(self.N - i -1)
        return res

    def compute(self, t):
        bernstein = self._berstein(t)
        curve_x = sum(self.polygon[i,0] * bernstein[i, :] \
                      for i in range(self.N))
        curve_y = sum(self.polygon[i,1] * bernstein[i, :] \
                      for i in range(self.N))
        return np.vstack((curve_x, curve_y,)).T


class FastBernsteinPolicy(BezierPolicy):
    import sys, os.path as path; sys.path.insert(0,path.dirname(path.realpath(__file__)))
    
    from binom import table

    del sys.path[0]

    @staticmethod
    def _fastexp(t,n):
        '''
        Create a table of (n+1)xlen(t) elements where:

            table[i] = t**i

        It makes extensive use of numpy views to improve performance
        '''

        r = np.empty( [n+1, len(t)] )
        r[0,:].fill(1); r[1,:] = t

        pi, i = 1, 1
        while i < n:
            newi = min(2*i+1,n)
            np.multiply(r[i,:], r[pi:pi + (newi-i)], # multiply this two
                        r[i+1:newi+1 ,:]) # output array
            pi=i; i=newi
        return r

    @classmethod
    def _binom(cls, n,k):
        '''
        binomial implementation using a table cache, preinitialised
        '''
        if n==-1 or k == -1: return 0
        if k > n-k: k = n-k
        index = ((n*(n+1)) >> 1) + k
        
        if len(cls.table) < index + k: 
            cls.table.extend( [-1]*( (((n+1)*(n+2)) >> 1) - len(cls.table))  )
        if cls.table[index] == -1: 
            cls.table[index] = cls._binom(n-1, k) + cls._binom(n-1, k-1)

        return cls.table[index]


    def _berstein(self,t):
        e = self._fastexp(t, self.N-1)

        res = np.empty((self.N, len(t),))
        for i in range(self.N):
            # carefully choose operations to not generate temporals
            np.multiply(e[i,:],e[self.N-i-1,::-1],res[i, :])
            res[i, :] *= self._binom(self.N - 1, i)
        return res


    def compute(self, t):
        bernstein = self._berstein(t)
        return np.dot(bernstein.T, self.polygon)
            

class FasterBernsteinPolicy(FastBernsteinPolicy):
    '''
    Fails because array index implies copy in numpy
    '''

    limit = 1 << 15
    precalc = FastBernsteinPolicy._fastexp(np.linspace(0,1,limit + 1), 1 << 9)
    
    @classmethod
    def _fastexp(cls,t,n):
        m = len(t) - 1 
        d, r = cls.limit / m , cls.limit % m
        def gen(end, step, off, points): 
            i,r,no = 0, 0, 0
            while r + no <= end: yield r + no; i += 1; r += step ; no = (i*off)/points
        return cls.precalc[:n+1,list(gen(cls.limit, d, r, m))]


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    poly = [ [0,0], [1,1.5],[2, 2], [1.25,0.7] ]

    ffbcurve = FasterBernsteinPolicy(poly)()
    fbcurve = FastBernsteinPolicy(poly)()
    bcurve = BernsteinPolicy(poly)()

    poly = np.array(poly)
    plt.figure(1)
    plt.subplot(311)
    plt.plot(bcurve[:,0], bcurve[:,1])
    plt.plot(poly[:,0], poly[:,1], "ro")

    plt.subplot(312)
    plt.plot(fbcurve[:,0], fbcurve[:,1])
    plt.plot(poly[:,0], poly[:,1], "ro")

    plt.subplot(313)
    plt.plot(ffbcurve[:,0], ffbcurve[:,1])
    plt.plot(poly[:,0], poly[:,1], "ro")


    plt.show()


    
    def eval_bezier(degree, t, algo):
        P = np.random.uniform(-20, 20, (degree + 1, 2))
        curve = algo(P)
        return curve #numpy array of size (num_points, 2)
    

    def eval_deCasteljau(degree, t):
        P = np.random.uniform(-20, 20, (degree + 1, 2))
    
        #enter here your computations
        curve = None #delete this line
        return curve #numpy array of size (num_points, 2)


    

    import timeit
    degree = 15
    num_points = 10000 
    number = 1000
    tt = np.linspace(0, 1, num_points)

    stuff = None

    def berns(P): BernsteinPolicy(P)(npoints = num_points)
    def fberns(P): FastBernsteinPolicy(P)(npoints = num_points)
    def ffberns(P): FasterBernsteinPolicy(P)(npoints = num_points)

    print(timeit.timeit("eval_bezier(degree, tt, ffberns)",
                        setup="from __main__ import eval_bezier, tt, degree, ffberns, FasterBernsteinPolicy", number=number))

    
    print(timeit.timeit("eval_bezier(degree, tt, fberns)",
                        setup="from __main__ import eval_bezier, tt, degree, fberns, FastBernsteinPolicy", number=number))


    print(timeit.timeit("eval_bezier(degree, tt, berns)",
                        setup="from __main__ import eval_bezier, tt, degree, berns, BernsteinPolicy", number=number))
    

    # print(timeit.timeit("eval_deCasteljau(degree, t)",
    #                     setup="from __main__ import eval_deCasteljau, degree, t, stuff"))
    

    
    
