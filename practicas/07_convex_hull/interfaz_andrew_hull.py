#!/usr/bin/env python2
# coding=utf-8

'''
Computational Geometry Assignments | (c) 2015 Pablo Cabeza & Diego González
license: [modified BSD](http://opensource.org/licenses/BSD-3-Clause)
'''

import numpy as np
import scipy.interpolate as sc
import matplotlib.pyplot as plt
from Tkinter import *
from ttk import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backend_bases import key_press_handler
from andrew_hull import convex_hull

class DrawPoints(object):
    def __init__(self, points):
        self.plt_points = points
        self.xs = list(self.plt_points.get_xdata())
        self.ys = list(self.plt_points.get_ydata())
        self.polygon = zip(self.xs, self.ys)
        self.plt_curve, = points.axes.plot([], [], 'k')
        self.selected_point = None

        self.cid_press = self.plt_points.figure.canvas.mpl_connect('button_press_event', self.press_event)
        self.cid_move = self.plt_points.figure.canvas.mpl_connect('motion_notify_event', self.move_event)
        self.cid_release = self.plt_points.figure.canvas.mpl_connect('button_release_event', self.release_event)
        self.cid_erease = self.plt_points.figure.canvas.mpl_connect('button_press_event', self.erease_event)

    def update_curve(self):
        p = convex_hull(np.array(self.polygon))
        self.plt_curve.set_data([x for x,_ in p], [y for _,y in p])
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
        if event.inaxes != self.plt_points.axes or event.button != 3:
            return
        state, indices = self.plt_points.contains(event)
        if state:
            index = indices['ind'][0]
            del self.xs[index]
            del self.ys[index]
            del self.polygon[index]
            self.update_curve()

def main(args=None):
    args = args or sys.argv

    #Se crea la ventana
    root = Tk()
    root.title("Andrew Hull")

    #Se crea la grafica
    fig = Figure(figsize=(5,4), dpi=100)
    ax = fig.add_subplot(111)
    ax.set_xlim([-10,10]), ax.set_ylim([-10,10])
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.show()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

    #Se generan los puntos e interpolamos
    points, = ax.plot([], [], 'ro')
    linebuilder = DrawPoints(points)

    #Se añade un pequeño tutorial:
    Ti = Text(root, height=3, width=50)
    Ti.pack()
    Ti.insert(END, '''INSTRUCTIONS: Left click to introduce points and  drag to move them or right click to erase them.   The convex hull is calculated autocatically.''')

    root.mainloop()

if __name__ == "__main__": sys.exit(main())
