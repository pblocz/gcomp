import timeit, sys, numpy as np


def eval_bezier(degree, t, algo):
    P = np.random.uniform(-20, 20, (degree + 1, 2))
    curve = algo(P)
    return curve #numpy array of size (num_points, 2)


class Tester(object):

    def __init__(self, fun, args, loops = 1000, extra_args = lambda: []):
        self.fun = fun
        self.loops = loops
        self.args = args
        self.extra_args = extra_args
    
    def time(self, name = ""):
        r = []
        for n in self.args:
            r.append(
                (n,
                timeit.timeit(stmt=lambda: self.fun(n, *self.extra_args()),
                              setup="pass",
                              number=self.loops),
                 )
            )
            print "%s with arg %s: %f" % (name, n, r[-1][1],) 
        return r

def bezier_tester(bezier):
    degree = 15
    def extra(): return [np.random.uniform(-20, 20, (degree + 1, 2))]
    return Tester(bezier, range(100,2000+1, 190), extra_args = extra)
                      
    
# Tester(fun_to_test, loops, extra_args = []).time()

def main(args = None):
    args = args or sys.argv
    
    from . import policies as plc

    t = bezier_tester(lambda n, p: plc.BernsteinPolicy(p)(n)).time("Bernstein")
    ft = bezier_tester(lambda n, p: plc.FastBernsteinPolicy(p)(n)).time("FastBerns")
    fft = bezier_tester(lambda n, p: plc.FasterBernsteinPolicy(p)(n)).time("FasterBerns") 

    import matplotlib.pyplot as plt

    fig = plt.figure()

    plt.plot([x for x,_ in t], [y for _,y in t], figure = fig)
    plt.plot([x for x,_ in ft], [y for _,y in ft], figure = fig)
    plt.plot([x for x,_ in fft], [y for _,y in fft], figure = fig)

    fig.savefig("algorithms.png")

    return (t,ft,fft)

if __name__ == "__main__": main()
