# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.special import binom
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D


class BezierPolicy(object):
    #Base class used to implement different curve methods
    def __init__(self, polygon):
        self.polygon = np.array(polygon)
        self.N = len(self.polygon)

    def __call__(self, npoints=100, *args, **kwargs): 
        return self.compute(np.linspace(0, 1, npoints), *args, **kwargs)
        

class BernsteinPolicy(BezierPolicy):

    def _berstein(self,t):
        res = np.zeros((self.N+1, len(t),))
        for i in range(self.N):        
            res[i,:] = binom(self.N-1, i)*t**i * (1-t)**(self.N-1-i)
        return res

    def compute(self, t):
        bernstein = self._berstein(t)
        curve_x = sum(self.polygon[i,0] * bernstein[i, :] for i in range(self.N))
        curve_y = sum(self.polygon[i,1] * bernstein[i, :] for i in range(self.N))
        return np.vstack((curve_x, curve_y,)).T


class DrawCurves(object):
    def __init__(self, points):
        self.plt_points = points
        self.xs = list(self.plt_points.get_xdata())
        self.ys = list(self.plt_points.get_ydata())
        self.polygon = zip(self.xs, self.ys)
        self.plt_curve, = points.axes.plot([], [])
        self.selected_point = None
        
        self.cid_press = self.plt_points.figure.canvas.mpl_connect('button_press_event', self.press_event)
        self.cid_move = self.plt_points.figure.canvas.mpl_connect('motion_notify_event', self.move_event)
        self.cid_release = self.plt_points.figure.canvas.mpl_connect('button_release_event', self.release_event)
        self.cid_erease = self.plt_points.figure.canvas.mpl_connect('button_press_event', self.erease_event)
    
    def update_Bezier(self):
        curve = BernsteinPolicy(self.polygon)(1000)
        self.plt_curve.set_data(curve[:,0], curve[:,1])
    
    def update_DeCasteljau(self):
        curve = BernsteinPolicy(self.polygon)(1000)
        self.plt_curve.set_data(curve[:,0], curve[:,1])

    def update_curve(self):
        #self.update_Bezier()
        self.update_DeCasteljau()
        self.plt_points.set_data(self.xs, self.ys)
        self.plt_points.figure.canvas.draw()
    
    def new_point(self, event):
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
        self.polygon.append([event.xdata,event.ydata])
        self.update_curve()
         
    def press_event(self, event):
        if event.inaxes != self.plt_points.axes or event.button != 1: 
            return
        state, indices = self.plt_points.contains(event)
        if not state:
            self.new_point(event)
        else:
            self.selected_point = indices['ind'][0]
            print self.selected_point
            
    def move_event(self, event):  
        if event.inaxes != self.plt_points.axes or self.selected_point is None:
            return
        self.xs[self.selected_point] = event.xdata
        self.ys[self.selected_point] = event.ydata
        self.polygon[self.selected_point] = [event.xdata,event.ydata]
        self.update_curve()     
        
    def release_event(self, event):  
        if event.inaxes != self.plt_points.axes: 
            return
        self.selected_point = None
        
    def erease_event(self, event):
        if event.inaxes != self.plt_points.axes or event.button != 3: 
            return
        state, indices = self.plt_points.contains(event)
        if state:
            index = indices['ind'][0]
            del self.xs[index]
            del self.ys[index]
            del self.polygon[index]
            self.update_curve() 
            
    



fig = plt.figure()
ax = fig.add_subplot(111)
points, = ax.plot([], [], 'ro')  # empty line
linebuilder = DrawCurves(points)

plt.show()