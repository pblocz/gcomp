#! /usr/bin/python -O
# -*- coding: utf-8 -*-
'''
Practica 3

@author: Pablo Cabeza García y Diego González
'''
import sys, re, argparse
import sympy as sp, numpy as np, matplotlib.pyplot as plt
from mayavi import mlab
from scipy.integrate import odeint
from scipy.spatial.distance import pdist

def hex2rgb(color): 
    return tuple(int(i,16)/255.0 for i in re.match('#?(.{2})(.{2})(.{2})',color).groups())

class BaseParametrized(object):
    def __init__(self, parameters):
        self.u, self.v = [v for v,_ in parameters]
        self.uR, self.vR = [list(r) for _,r in parameters]

    @property 
    def parameters(self): return [(self.u, self.uR), (self.v, self.vR)]
    
    @property 
    def variables(self): return [self.u, self.v]


    @staticmethod 
    def _range(R,points): return np.linspace(*(R + [points]))

    @property 
    def urange(self, npoints = 100): return self._range(self.uR, npoints)

    @property 
    def vrange(self, npoints = 100): return self._range(self.vR, npoints)



class Geodesic(object):
    def __init__(self, points):
        self.points = points
        self.figure = None

    def plot(self, figure=None):
        self.figure = figure or self.figure or plt.figure()
        plt.plot(self.points[:,0],self.points[:,1], figure = self.figure)
        return self.figure

    def show(self): plt.show(); return self

    def clear(self): 
        self.figure.clear()
        self.figure = None
        return self

    def savefig(self,*args,**kwargs): 
        plt.savefig(*args,figure= self.figure,**kwargs)
        return self


class FFForm(BaseParametrized):
    def __init__(self,E,F,G, parameters=None):
        super(FFForm, self).__init__(parameters)
        self.E, self.F, self.G = E, F, G

    def _curvize(self,surf, fvs): return surf.subs(zip(self.variables,fvs))

    def _rhs_gen(self):
        u, v, t = self.u, self.v, sp.symbols('t')
        du,dv = sp.symbols('du dv')
        U = sp.Matrix((sp.Function('u')(t),sp.Function('v')(t),))
        U_ = U.diff(t)

        I = sp.Matrix([[self.E,self.F],[self.F,self.G]])
        Iu = self._curvize(I.diff(u), list(U))
        Iv = self._curvize(I.diff(v),  list(U))
        I = self._curvize(I, list(U))

        U__ = (1/2.0 * sp.Matrix([U_.T * Iu * U_, U_.T*Iv*U_]).T - U_.T*(Iu*U_[0] + Iv*U_[1]))\
              * (I**-1); 
        U__ = sp.simplify(U__).T

        MU__ = U__.subs(zip(list(U_), [du,dv] ,))
        MU__ = MU__.subs(zip(list(U), [u,v],))

        return sp.lambdify((u,v,du,dv), sp.Matrix([du,dv,MU__[0],MU__[1]]).T,
                           modules='numpy')


    def geodesic(self, p0, vt, interval, npoints=100):
        rhs = self._rhs_gen()
        def rhs_eqs(Y,t): return np.squeeze(np.asarray(rhs(*Y)))
        init_cond = list(p0) + list(vt)

        I0,In = interval
        solv = odeint(rhs_eqs, init_cond, np.linspace(I0,In,npoints))[:,0:2]
        return Geodesic(solv)


class Parametrization3D(BaseParametrized):
    def __init__(self, x,y,z,parameters=None):
        super(Parametrization3D,self).__init__(parameters)

        self.x, self.y, self.z = x, y, z
        self.lx, self.ly, self.lz = [sp.lambdify(self.variables,s,modules='numpy') \
                                     for s in [x,y,z]]
        self.figure = {}

    @property 
    def surface(self): return sp.Matrix([self.x, self.y, self.z])

    @property 
    def lsurface(self): return [self.lx, self.ly, self.lz]


    def firstForm(self):
        s_u, s_v = self.surface.diff(self.u), self.surface.diff(self.v)
        E,F,G = s_u.dot(s_u.T), s_u.dot(s_v.T), s_v.dot(s_v.T)
        return FFForm(E,F,G, self.parameters)


    def _mesh(self, points,**kwargs):
        u_mesh, v_mesh = np.meshgrid(self.urange, self.urange)

        mlab.mesh(*[k(u_mesh,v_mesh) for k in self.lsurface],figure = self.figure,
                  **kwargs)
    
    def _plot3d(self,curve, points,**kwargs):
        mlab.plot3d(*[(k(curve[:,0],curve[:,1])) for k in self.lsurface], figure= self.figure,
                    **kwargs)

    def plot(self,curve = None, points = 100, figure=None,**kwargs):
        self.figure = figure or self.figure or mlab.figure()

        if curve: self._plot3d(curve.points,points, **kwargs)
        else: self._mesh(points, **kwargs)

        return self.figure

    def show(self,figure=None): mlab.draw(figure = figure or self.figure); mlab.show(); return self
    def clear(self): 
        try: mlab.clf(self.figure); mlab.close(self.figure)
        except: pass

        self.figure = None
        return self
    def savefig(self,*args,**kwargs): mlab.savefig(*args,figure=self.figure,**kwargs); return self


def modular_plot(geo):
    "plots the geodesic in the space [0,2pi]x[0,2pi]"

    fig = plt.figure()
    p = geo.points; pi2 = 2*np.pi
    ul = p[:,0]; l = []

    Mul, mul = np.amax(ul), np.amin(ul)
    Mupi, mupi = np.ceil(Mul / (2*np.pi)), np.floor(mul / (2*np.pi))
    for i in np.arange(pi2*mupi,pi2*Mupi,pi2):
        ull = p[(i < p[:,0]) & (p[:,0] < i+pi2)]

        vll = ull[:,1]; Mvl, mvl = np.amax(vll), np.amin(vll)
        Mvpi, mvpi = np.ceil(Mvl / pi2), np.floor(mvl / pi2)
        for j in np.arange(pi2*mvpi,pi2*Mvpi,pi2):
            vec = ull[(j < ull[:,1]) & ( ull[:,1] < j+pi2)] % pi2
            if vec.size != 0: l.append(vec)

    for d in l: plt.plot(d[:,0],d[:,1], figure=fig)
    return fig



def main():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-show','-n', dest='show', action="store_false", default=True)
    args = parser.parse_args()

    # Define variables to use
    u,v = sp.symbols('u,v')    

    # Define torus and get first fundamental form
    R = 2; r = 1
    torus = Parametrization3D(sp.cos(u) * (R + r*sp.cos(v)),
                              (sp.sin(u) * (R + r*sp.cos(v))),
                              r *sp.sin(v),
                              parameters = [(u,(0,2*np.pi)),(v,(0,2*np.pi))])
    torusff = torus.firstForm()

    
    # Example for generators circunferences #
    torus.plot()
    torus.plot(torusff.geodesic((0,0),(0,1), (0,2*np.pi)))
    torus.plot(torusff.geodesic((0,0),(1,0), (0,2*np.pi)), color=hex2rgb('#A753A4'))
    torus.savefig('generators.png', size=(1000,1000))
    if args.show: torus.show()

    # Example for perodic geodesic #
    periodic = torusff.geodesic([0,np.pi/4.0],[1,1],[0,100*np.pi],npoints=10000)
    torus.clear().plot()
    torus.plot(periodic, line_width=1.0, tube_radius=None) # this disables big tubes
    torus.savefig('periodic.png', size=(1000,1000))
    if args.show: torus.show()

    fig = modular_plot(periodic)
    fig.savefig('periodic_uv.png')
    if args.show: fig.show()
    fig.clear()


    # Example for non periodic geodesic #
    minidense = torusff.geodesic([0,3*np.pi/4.0],[1,np.sqrt(2)],[0,500*np.pi],npoints=100000)
    torus.clear().plot()
    torus.plot(minidense, line_width=1.0, tube_radius=None) # this disables big tubes
    torus.savefig('minidense.png', size=(1000,1000))
    if args.show: torus.show()

    fig = modular_plot(minidense)
    fig.savefig('minidense_uv.png')
    if args.show: fig.show()
    fig.clear()

    # Same Example, but more loops #
    dense = torusff.geodesic([0,3*np.pi/4.0],[1,np.sqrt(2)],[0,5000*np.pi],npoints=1000000)
    torus.clear().plot()
    torus.plot(dense, line_width=1.0, tube_radius=None) # this disables big tubes
    torus.savefig('dense.png', size=(1000,1000))
    if args.show: torus.show()
    
    fig = modular_plot(dense)
    fig.savefig('dense_uv.png')
    if args.show: fig.show()
    fig.clear()


    # Poincare half-plane geodesics #
    poincare = FFForm(1/v**2, 0, 1/v**2, parameters = [(u, []), (v, [])])
    
    fig = None
    for exp in [i*2 for i in range(1,11)]:
        fig = poincare.geodesic((0,1),(1,exp),(-100,100),npoints=10000)\
                      .plot(figure = fig)
        fig = poincare.geodesic((0,1),(-1,exp),(-100,100),npoints=10000)\
                      .plot(figure = fig)

    fig = poincare.geodesic((0,40),(1,0),(0,1000),npoints=10000)\
                  .plot(figure = fig)
    fig = poincare.geodesic((0,40),(-1,0),(0,1000),npoints=10000)\
                  .plot(figure = fig)

    plt.axes().set_aspect('equal')
    plt.savefig('poincare.png', figure=fig)
    if args.show: plt.show()

    return 0

if __name__ == "__main__": sys.exit(main())
