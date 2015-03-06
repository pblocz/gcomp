# -*- coding: utf-8 -*-
'''
Practica 2

@author: Pablo Cabeza García y Diego González
'''
import sympy as sp, numpy as np
from mayavi import mlab
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist

def geodesic(surface, p0, vt, Interv, t0, npoints = 300, varis=None):
    '''
    surface: (u,v) -> R3
    p0: (u0,v0) | x0 = surface(u0,v0)
    vt: xu(p0)u0' + xv(p0)v0'
    I: [a,b]
    t0: t0 in I
    '''
    varis = varis or surface.free_symbols
    u,v = varis
    du,dv = sp.symbols('du dv')
    t = sp.symbols('t')
    U = sp.Matrix((sp.Function('u')(t),sp.Function('v')(t),))
    U_ = U.diff(t)

    def curvize(surf,fun=list(U),varis = varis): return surf.subs(zip(varis, fun))

    s_u, s_v = surface.diff(u), surface.diff(v)
    E,F,G = s_u.dot(s_u.T), s_u.dot(s_v.T), s_v.dot(s_v.T)

    I = sp.Matrix([[E,F],[F,G]])
    Iu = curvize(I.diff(u))
    Iv = curvize(I.diff(v))
    I = curvize(I)

    
    U__ = (1/2.0 * sp.Matrix([U_.T * Iu * U_, U_.T*Iv*U_]).T - U_.T*(Iu*U_[0] + Iv*U_[1])) * (I**-1); U__ = sp.simplify(U__.T)

    U__ = U__.T
    # sp.pprint(U__)

    MU__ = U__.subs(zip(list(U_), [du,dv] ,))
    MU__ = MU__.subs(zip(list(U), [u,v],))
    # sp.pprint(MU__)

    ## NUMPY ##
    # Resolver la ode
    rhs = sp.lambdify((u,v,du,dv),sp.Matrix([du,dv,MU__[0],MU__[1]]).T, modules='numpy')
    def rhs_eqs(Y,t): return np.squeeze(np.asarray(rhs(*Y)))
    init_cond = p0 + vt

    I0,In = Interv
    # intervs = [np.linspace(t0,I0,npoints), np.linspace(t0,In, npoints)[1:]]
    # solvs = [odeint(rhs_eqs, init_cond, interv) for interv in intervs]
    # solv = np.concatenate((solvs[0][::-1],solvs[1]),axis=0)


    solv = odeint(rhs_eqs, init_cond, np.linspace(I0,In,npoints))
    
    x = sp.lambdify((u,v),surface[0],modules='numpy')
    y = sp.lambdify((u,v),surface[1],modules='numpy')
    z = sp.lambdify((u,v),surface[2],modules='numpy')
    
    curve_u = solv[:,0] % 2*np.pi; curve_v = solv[:,1] % 2*np.pi

    # print zip(curve_u, curve_v)
    # plt.axes().set_aspect('equal')
    # plt.plot(curve_u, curve_v)
    # plt.show()


    print pdist([(curve_u[0],curve_v[0]), (curve_u[-1],curve_v[-1])])

    u_range = np.linspace(0,2*np.pi,100)
    v_range = np.linspace(0,2*np.pi,100)

    u_mesh, v_mesh = np.meshgrid(u_range, v_range)

    def unique(a):
        a = np.sort(a)
        b = np.diff(a)
        b = np.r_[1, b]
        return a[b != 0]


    plt.plot(curve_u,curve_v)
    plt.show()

    # print [k(curve_u,curve_v) for k in [x, y, z]]
    mlab.mesh(*[k(u_mesh,v_mesh) for k in [x, y, z]])
    mlab.plot3d(*[(k(curve_u,curve_v)) for k in [x, y, z]])
    mlab.show()

    
u,v = sp.symbols('u,v')
# surf = sp.Matrix((sp.cosh(u) * sp.cos(v),
#         sp.cosh(u) * sp.sin(v),
#         sp.sinh(u),))


# geodesic(surf,[0,0],[0,1],[np.pi,2*np.pi-0.1],np.pi)

R = 2; r = 1
surf = sp.Matrix((sp.cos(u) * (R + r*sp.cos(v)),
                 (sp.sin(u) * (R + r*sp.cos(v))),
                  r *sp.sin(v),))


# Circunferencia circulo generador (grosor del donuts)
# geodesic(surf,[0,0],[0,1],[0,2*np.pi],1,varis=[u,v])

# # Circunferencia al rededor del círculo de revolución (radio del donuts)
# geodesic(surf,[0,0],[1,0],[0,2*np.pi],1,varis=[u,v])

# Geodésica periódica con muchos puntos
# geodesic(surf,[0,0],[1,np.pi],[0,2*np.pi],1,varis=[u,v],npoints=10000)
# geodesic(surf,[0,np.pi/4.0],[1,1],[0,100*np.pi],1,varis=[u,v],npoints=10000)

# geodesic(surf,[0,0],[1,4],[0,40*np.pi],1,varis=[u,v],npoints=10000)
# geodesic(surf,[0,0],[1,1000],[0,5*np.pi],1,varis=[u,v],npoints=1000)

geodesic(surf,[0,np.pi/4.0],[1,np.sqrt(2)],[0,1000*np.pi],1,varis=[u,v],npoints=1000000)

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
