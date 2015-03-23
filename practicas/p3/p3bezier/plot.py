# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.special import binom
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D

from p3bezier import policies as plc


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
        curve = plc.FastBernsteinPolicy(self.polygon)(1000)
        self.plt_curve.set_data(curve[:,0], curve[:,1])
    
    def update_DeCasteljau(self):
        curve = plc.DeCasteljauFastPolicy(self.polygon)(1000)
        self.plt_curve.set_data(curve[:,0], curve[:,1])

    def update_curve(self):
        "update the drawn curve"
        self.update_Bezier()
        # self.update_DeCasteljau()
        self.plt_points.set_data(self.xs, self.ys)
        self.plt_points.figure.canvas.draw()
    
    def new_point(self, event):
        "adds a point when the canvas is clicked"
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
        self.polygon.append([event.xdata,event.ydata])
        self.update_curve()
         
    def press_event(self, event):
        "Event when press with primary key (add or move point)"
        if event.inaxes != self.plt_points.axes or event.button != 1: 
            return
        state, indices = self.plt_points.contains(event)
        if not state:
            self.new_point(event)
        else:
            self.selected_point = indices['ind'][0]
            # print self.selected_point
            
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
        "Event when pressed with secondary key (delete point)"
        if event.inaxes != self.plt_points.axes or event.button != 3: 
            return
        state, indices = self.plt_points.contains(event)
        if state:
            index = indices['ind'][0]
            del self.xs[index]
            del self.ys[index]
            del self.polygon[index]
            self.update_curve() 
            
def main():
    print "click and move some points to see the bezier curve generated"
    print "use secondary click to delete points"

    fig = plt.figure()
    ax = fig.add_subplot(111)
    points, = ax.plot([], [], 'ro')  # empty line
    linebuilder = DrawCurves(points)

    plt.show()

if __name__ == '__main__': main()
