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

    
class DrawPoints(object):
    def __init__(self, points):
        self.plt_points = points
        self.xs = list(self.plt_points.get_xdata())
        self.ys = list(self.plt_points.get_ydata())
        self.polygon = zip(self.xs, self.ys)
        self.plt_curve, = points.axes.plot([], [])
        self.selected_point = None
        
        self.knots = 'otro'
        self.method = 'newton'
        self.L = 0
        self.libraries = False
        self.num_points = 100
        
        self.cid_press = self.plt_points.figure.canvas.mpl_connect('button_press_event', self.press_event)
        self.cid_move = self.plt_points.figure.canvas.mpl_connect('motion_notify_event', self.move_event)
        self.cid_release = self.plt_points.figure.canvas.mpl_connect('button_release_event', self.release_event)
        self.cid_erease = self.plt_points.figure.canvas.mpl_connect('button_press_event', self.erease_event)
    
    def update_data(knots, method, L, libraries, num_points):
        self.knots = knots
        self.method = method
        self.L = L
        self.libraries = libraries
        self.num_points = num_points
        
    
    def update_curve(self):
        #self.updatePolynom()
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
    root.title("Intersección de curvas de Bezier")
    
    #Se crea la zona para seleccionar los nodos
    knots = Frame(root)
    knots.pack(anchor=W) 
    klabel = Label(knots, text="KNOTS:                ")
    klabel.pack(side=LEFT, padx=2, pady=2)
    vn = BooleanVar()
    vn.set(False)
    #N = StringVar()
    #N.set("10")
    R1=Radiobutton(knots, text="Chebychev   ", variable=vn, value=True)
    R2=Radiobutton(knots, text="Otro    ", variable=vn, value=False)
    #Nvalue = Entry(knots, textvariable=N, width=3, justify=CENTER)
    R1.pack(anchor=W, side=LEFT)
    R2.pack(anchor=W, side=LEFT)
    #Nvalue.pack()
    
    #Se crea la zona para seleccionar el método
    method = Frame(root)
    method.pack(anchor=W) 
    mlabel = Label(method, text="METHOD:            ")
    mlabel.pack(side=LEFT, padx=2, pady=2)
    vm = BooleanVar()
    vm.set(True)
    L = StringVar()
    L.set("0")
    Lvalue = Entry(method, textvariable=L, width=3, justify=CENTER, state='disabled')
    def check():
        if vm.get() == 0:
            Lvalue.configure(state='normal')
        else:
            Lvalue.configure(state='disabled')
    R3=Radiobutton(method, text="Polinomio de Newton   ", variable=vm, value=True, command=check)
    R4=Radiobutton(method, text="Minimos Cuadrados,    L=", variable=vm, value=False, command=check)
    R3.pack(anchor=W, side=LEFT)
    R4.pack(anchor=W, side=LEFT)
    Lvalue.pack()
    
    #Se crea la zona para seleccionar si se usan librerias o no
    library = Frame(root)
    library.pack(anchor=W) 
    llabel = Label(library, text="LIBRARY:              ")
    llabel.pack(side=LEFT, padx=2, pady=2)
    vl = BooleanVar()
    vl.set(False)
    R5=Radiobutton(library, text="SIN librerias   ", variable=vl, value=False, command=check)
    R6=Radiobutton(library, text="CON librerias", variable=vl, value=True, command=check)
    R5.pack(anchor=W, side=LEFT)
    R6.pack(anchor=W, side=LEFT)
    
    #Se crea la zona para seleccionar ver el numero de puntos:
    points = Frame(root)
    points.pack(anchor=W) 
    plabel = Label(points, text="NUM. POINTS:     ")
    plabel.pack(side=LEFT, padx=2, pady=2)
    P = StringVar()
    P.set("100")
    Pvalue = Entry(points, textvariable=P, width=5, justify=CENTER)
    Pvalue.pack()
        
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
    def function():
        if vn.get() == 1:
            knots = 'chebyshev'
        else:
            knots = 'otro'
        #nk = int(Nvalue.get())
        if vm.get():
            method = 'least_squares'
            L = int(Lvalue.get())
        else:
            method = 'newton'
            L = 0
        libraries = vl.get()
        num_points = int(Pvalue.get())
        
        linebuilder.update_data(knots, method, L, libraries, num_points)
        linebuilder.update_curve()
                
        print "INTERPOLADO"
    button = Button(root, text="INTERPOLAR", command=function)
    button.pack()
        
    root.mainloop()

if __name__ == "__main__": sys.exit(main())