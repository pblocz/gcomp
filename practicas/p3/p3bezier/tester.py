import timeit, sys, numpy as np
import matplotlib.pyplot as plt

def eval_bezier(degree, t, algo):
    P = np.random.uniform(-20, 20, (degree + 1, 2))
    curve = algo(P)
    return curve #numpy array of size (num_points, 2)


class Tester(object):
    '''Tester class for bezier curves'''

    def __init__(self, fun, args, loops = 1000, extra_args = lambda: []):
        self.fun = fun
        self.loops = loops
        self.args = args
        self.extra_args = extra_args
    
    def time(self, name = ""):
        r = []
        for n in self.args:
            t = timeit.timeit(stmt=lambda: self.fun(n, *self.extra_args()),
                              setup="pass",
                              number=self.loops)
            r.append((n, t, t / float(self.loops)))
            print "%s with arg %s: %f -> loop %f" % (name,n, r[-1][1], r[-1][2]) 
        self.last_time = r
        self.last_name = name
        return r

    def loop_time(self, *args , **kwargs):
        r = self.time(*args , **kwargs)
        self.last_time = [(n,k/float(self.loops),) for n,k, in r]

        return self.last_time

    def plot(self, **kwargs):
        return plt.plot([x for x,_,_ in self.last_time], 
                        [y for _,y,_ in self.last_time], 
                        label=self.last_name, **kwargs)

    def plot_loops(self, **kwargs):
        return plt.plot([x for x,_,_ in self.last_time], 
                        [y for _,_,y in self.last_time], 
                        label=self.last_name, **kwargs)



def bezier_tester(bezier,args = range(100,2000+1, 190), loops = 1000):
    degree = 15
    def extra(): return [np.random.uniform(-20, 20, (degree + 1, 2))]
    return Tester(bezier, args, extra_args = extra, loops = loops)
                      
    

def main(args = None):
    args = args or sys.argv
    
    from . import policies as plc
    from p3cbezier import bernstein as CBernstein

    print "Basic test (as published in the campus)"

    ct, = bezier_tester(lambda n, p: CBernstein(p,n),[100], loops = 10000).time("CBerns") 
    t, = bezier_tester(lambda n, p: plc.BernsteinPolicy(p)(n),[100], loops = 10000).time("Bernstein")
    ft, = bezier_tester(lambda n, p: plc.FastBernsteinPolicy(p)(n),[100], loops = 10000).time("FastBerns")
    fft, = bezier_tester(lambda n, p: plc.FasterBernsteinPolicy(p)(n),[100], loops = 10000).time("FasterBerns") 
    djt, = bezier_tester(lambda n, p: plc.DeCasteljauPolicy(p)(n),[100], loops = 10000).time("DeCasteljau") 
    djtf, = bezier_tester(lambda n, p: plc.DeCasteljauFastPolicy(p)(n),[100], loops = 10000).time("DeCasteljauFast")

    print
    print "Press any key to continue to a more complex test, ctrl+c to quit"
    raw_input()

    # ct = bezier_tester(lambda n, p: CBernstein(p,n)).time("CBerns") 
    t = bezier_tester(lambda n, p: plc.BernsteinPolicy(p)(n)).time("Bernstein")
    ft = bezier_tester(lambda n, p: plc.FastBernsteinPolicy(p)(n)).time("FastBerns")
    fft = bezier_tester(lambda n, p: plc.FasterBernsteinPolicy(p)(n)).time("FasterBerns") 
    djt = bezier_tester(lambda n, p: plc.DeCasteljauFastPolicy(p)(n)).time("DeCasteljau") 


    fig = plt.figure()

    # cl, = plt.plot([x for x,_ in ct], [y for _,y in ct], figure = fig, label="cbernstein")
    l, = plt.plot([x for x,_ in t], [y for _,y in t], figure = fig, label="Bernstein")
    fl, = plt.plot([x for x,_ in ft], [y for _,y in ft], figure = fig, label="FastBernstein")
    ffl, = plt.plot([x for x,_ in fft], [y for _,y in fft], figure = fig, label="FasterBernstein")
    dj, = plt.plot([x for x,_ in djt], [y for _,y in djt], figure = fig, label="De Casteljau")

    plt.legend([l,fl,ffl,dj,])
    fig.savefig("algorithms-pablodiego.png")

    print "Saved images in 'algorithms-pablodiego.png'"
    plt.show()

    # return (t, ft, fft,dj)

if __name__ == "__main__": main()
