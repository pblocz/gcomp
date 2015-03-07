# -*- coding: utf-8 -*-
'''
Practica 3

@author: Pablo Cabeza García y Diego González
'''
import sys
import sympy as sp, numpy as np, matplotlib.pyplot as plt
from mayavi import mlab
from scipy.integrate import odeint
from scipy.spatial.distance import pdist

class BaseParametrized(object): pass


class Geodesic(object):
    def __init__(self, points):
        self.points = points

    def plot(self):
        print self.points
        return plt.plot(self.points[:,0],self.points[:,1])
        #return plt.plot(self.points.T)
    def show(self): 
        plt.show()


class FFForm(object):
    def __init__(self,E,F,G, parameters=None):
        self.E, self.F, self.G = E, F, G

        variables = [v for v,_ in parameters]
        self.u, self.v = variables
        self.uR, self.vR = [list(r) for _,r in parameters]

    def _curvize(self,surf, fvs): return surf.subs(zip([self.u,self.v],fvs))

    def _rhs_gen(self):
        u, v, t = self.u, self.v, sp.symbols('t')
        du,dv = sp.symbols('du dv')
        U = sp.Matrix((sp.Function('u')(t),sp.Function('v')(t),))
        U_ = U.diff(t)

        I = sp.Matrix([[self.E,self.F],[self.F,self.G]])
        Iu = self._curvize(I.diff(u), list(U))
        Iv = self._curvize(I.diff(v),  list(U))
        I = self._curvize(I, list(U))

        U__ = (1/2.0 * sp.Matrix([U_.T * Iu * U_, U_.T*Iv*U_]).T - U_.T*(Iu*U_[0] + Iv*U_[1])) * (I**-1); 
        U__ = sp.simplify(U__).T

        MU__ = U__.subs(zip(list(U_), [du,dv] ,))
        MU__ = MU__.subs(zip(list(U), [u,v],))

        return sp.lambdify((u,v,du,dv),sp.Matrix([du,dv,MU__[0],MU__[1]]).T, modules='numpy')


    def geodesic(self, p0, vt, interval, npoints=100):
        rhs = self._rhs_gen()
        def rhs_eqs(Y,t): return np.squeeze(np.asarray(rhs(*Y)))
        init_cond = list(p0) + list(vt)

        I0,In = interval
        solv = odeint(rhs_eqs, init_cond, np.linspace(I0,In,npoints))[:,0:2]
        return Geodesic(solv)


class Parametrization3D(object):

    def __init__(self, x,y,z,parameters=None):
        variables = [v for v,_ in parameters]
        self.u, self.v = variables
        self.uR, self.vR = [list(r) for _,r in parameters]

        self.x, self.y, self.z = x, y, z
        self.lx, self.ly, self.lz = [sp.lambdify(variables,s,modules='numpy') for s in [x,y,z]]
        self.figure = {}


    @property
    def parameters(self): return [(self.u, self.uR,), (self.v, self.vR,)]

    @property
    def surface(self): return sp.Matrix([self.x, self.y, self.z])

    @property
    def lsurface(self): return [self.lx, self.ly, self.lz]


    def firstForm(self):
        s_u, s_v = self.surface.diff(self.u), self.surface.diff(self.v)
        E,F,G = s_u.dot(s_u.T), s_u.dot(s_v.T), s_v.dot(s_v.T)
        return FFForm(E,F,G, self.parameters)


    def _mesh(self, points):
        u_range = np.linspace(*(self.uR + [points]))
        v_range = np.linspace(*(self.uR + [points]))
        u_mesh, v_mesh = np.meshgrid(u_range, v_range)

        mlab.mesh(*[k(u_mesh,v_mesh) for k in self.lsurface],figure = self.figure)
    
    def _plot3d(self,curve, points):
        u_range = np.linspace(*(self.uR + [points]))
        v_range = np.linspace(*(self.uR + [points]))
        mlab.plot3d(*[(k(curve[:,0],curve[:,1])) for k in self.lsurface], figure= self.figure)

    def plot(self,curve = None, points = 100, figure=None):
        self.figure = figure or self.figure or mlab.figure()

        if curve: self._plot3d(curve.points,points)
        else: self._mesh(points)

        return self.figure

    def show(self,figure=None): mlab.draw(figure = figure or self.figure); mlab.show()

    
def main():
    
    u,v = sp.symbols('u,v')
    
    R = 2; r = 1
    para = Parametrization3D( sp.cos(u) * (R + r*sp.cos(v)),
                              (sp.sin(u) * (R + r*sp.cos(v))),
                              r *sp.sin(v),
                              parameters = [(u,(0,2*np.pi)),(v,(0,2*np.pi))])
    para.plot()
    geo = para.firstForm().geodesic((0,0),(0,1),(0,2*np.pi))
    para.plot(geo)
    para.show()

    geo.plot()
    geo.show()


    poincare = FFForm(1/v**2, 0, 1/v**2, parameters = [(u, []), (v, [])])
    pgeo = poincare.geodesic((1,10),(1,10),(-100,100),npoints=10000)
    pgeo.plot()
    pgeo.show()

    return 0

if __name__ == "__main__": sys.exit(main())

# Circunferencia circulo generador (grosor del donuts)
# geodesic(surf,[0,0],[0,1],[0,2*np.pi],1,varis=[u,v])

# # Circunferencia al rededor del círculo de revolución (radio del donuts)
# geodesic(surf,[0,0],[1,0],[0,2*np.pi],1,varis=[u,v])

# Geodésica periódica con muchos puntos
# geodesic(surf,[0,0],[1,np.pi],[0,2*np.pi],1,varis=[u,v],npoints=10000)
# geodesic(surf,[0,np.pi/4.0],[1,1],[0,100*np.pi],1,varis=[u,v],npoints=10000)

# geodesic(surf,[0,0],[1,4],[0,40*np.pi],1,varis=[u,v],npoints=10000)
# geodesic(surf,[0,0],[1,1000],[0,5*np.pi],1,varis=[u,v],npoints=1000)

# geodesic(surf,[0,np.pi/4.0],[1,np.sqrt(2)],[0,1000*np.pi],1,varis=[u,v],npoints=1000000)

# theta = np.linspace(0,2*np.pi,100)
# t = np.linspace(-1,1,100)

# theta_mesh, t_mesh = np.meshgrid(theta,t)

# x = np.cosh(t_mesh) * np.cos(theta_mesh)
# y = np.cosh(t_mesh) * np.sin(theta_mesh)
# z = np.sinh(t_mesh)

# x_ = np.cosh(t) * np.cos(theta)
# y_ = np.cosh(t) * np.sin(theta)
# z_ = np.sinh(t)



# mlab.plot3d(x_, y_, z_)
# mlab.mesh(x,y,z)
# mlab.show()
