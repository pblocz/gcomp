#! /usr/bin/env python2
# -*- encoding: utf-8 -*-

'''
authors: Pablo Cabeza & Diego González
'''

import numpy as np
import scipy.interpolate as sc
import matplotlib.pyplot as plt
from Tkinter import *
from ttk import * 
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backend_bases import key_press_handler
from splines import spline2d
  
class DrawPoints(object):
    def __init__(self, points):
        self.plt_points = points
        self.xs = list(self.plt_points.get_xdata())
        self.ys = list(self.plt_points.get_ydata())
        self.polygon = zip(self.xs, self.ys)
        self.plt_curve, = points.axes.plot([], [], 'k')
        self.selected_point = None
        
        self.a = 0
        self.b = 1
        self.xi = []
        self.k = 1
        self.nu = []
        self.A = []
        self.num_dots = 100
        
        self.cid_press = self.plt_points.figure.canvas.mpl_connect('button_press_event', self.press_event)
        self.cid_move = self.plt_points.figure.canvas.mpl_connect('motion_notify_event', self.move_event)
        self.cid_release = self.plt_points.figure.canvas.mpl_connect('button_release_event', self.release_event)
        self.cid_erease = self.plt_points.figure.canvas.mpl_connect('button_press_event', self.erease_event)
        
    def update_data(self, a, b, xi, k, nu, num_dots):
        self.a = a
        self.b = b
        self.xi = xi
        self.k = k
        self.nu = nu
        self.A = self.plt_points
        self.num_dots = 100
              
    def update_curve(self):
        p = spline2d(self.a, self.b, self.xi, self.k, self.nu, self.A, self.num_dots)
        self.plt_curve.set_data(p[:,0],p[:,1])
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
    root.title("Plane Spline Curve")
    
    #Se crea la zona para seleccionar los nodos
    knots = Frame(root)
    knots.pack(anchor=W) 
    alabel = Label(knots, text="a: ")
    alabel.pack(side=LEFT, padx=2, pady=2)
    La = StringVar()
    La.set("0")
    Lavalue = Entry(knots, textvariable=La, width=3, justify=CENTER)
    Lavalue.pack(side=LEFT)
    blabel = Label(knots, text="    b: ")
    blabel.pack(side=LEFT, padx=2, pady=2)
    Lb = StringVar()
    Lb.set("1")
    Lbvalue = Entry(knots, textvariable=Lb, width=3, justify=CENTER)
    Lbvalue.pack(side=LEFT)
    xilabel = Label(knots, text="    List of breakpoints: ")
    xilabel.pack(side=LEFT, padx=2, pady=2)
    Lxi = StringVar()
    Lxi.set("")
    Lxivalue = Entry(knots, textvariable=Lxi, width=30, justify=LEFT)
    Lxivalue.pack(side=LEFT)
    
    #Se crea la zona para seleccionar el orden
    order = Frame(root)
    order.pack(anchor=W) 
    klabel = Label(order, text="Order of the curve: ")
    klabel.pack(side=LEFT, padx=2, pady=2)
    Lk = StringVar()
    Lk.set("2")
    Lkvalue = Entry(order, textvariable=Lk, width=3, justify=CENTER)
    Lkvalue.pack(side=LEFT)
    
    #Se crea la zona para seleccionar el numero de condiciones
    smooth = Frame(root)
    smooth.pack(anchor=W) 
    nulabel = Label(smooth, text="Number of smoothness conditions: ")
    nulabel.pack(side=LEFT, padx=2, pady=2)
    Lnu = StringVar()
    Lnu.set("")
    Lnuvalue = Entry(smooth, textvariable=Lnu, width=30, justify=LEFT)
    Lnuvalue.pack(side=LEFT)
    
    #Se crea la zona para seleccionar el numero de puntos:
    points = Frame(root)
    points.pack(anchor=W) 
    plabel = Label(points, text="Number of dots: ")
    plabel.pack(side=LEFT, padx=2, pady=2)
    P = StringVar()
    P.set("100")
    Pvalue = Entry(points, textvariable=P, width=5, justify=CENTER)
    Pvalue.pack()
    buttons = Frame(root)
    buttons.pack() 
            
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
    points, = ax.plot([], [], 'ro-')
    linebuilder = DrawPoints(points)
    
    #Se crea el boton para interpolar
    def update():
        a = float(Lavalue.get())
        b = float(Lbvalue.get())
        xi = np.fromstring(Lxivalue.get(),dtype=float, sep=',')
        k = int(Lkvalue.get())
        nu = np.fromstring(Lnuvalue.get(),dtype=float, sep=',')
        num_dots = int(Pvalue.get())
        linebuilder.update_data(a, b, xi, k, nu, num_dots)
        print a, b, xi, k, nu, num_dots
        #linebuilder.update_curve()
    button = Button(buttons, text="UPDATE", command=update)
    button.pack()
    
    #Se añade un pequeño tutorial:
    Ti = Text(root, height=4, width=50)
    Ti.pack()
    Ti.insert(END, '''INSTRUCTIONS:  Left click to introduce points and find the plane spline curve. You can drag points  to move them or right click to erase them. Configure above options to see different curves.''')
    
    root.mainloop()

if __name__ == "__main__": sys.exit(main())