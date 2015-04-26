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
from polynomial_curve_fitting import polynomial_curve_fitting
  
class DrawPoints(object):
    def __init__(self, points):
        self.plt_points = points
        self.xs = list(self.plt_points.get_xdata())
        self.ys = list(self.plt_points.get_ydata())
        self.polygon = zip(self.xs, self.ys)
        self.plt_curve, = points.axes.plot([], [], 'k')
        self.selected_point = None
        
        self.knots = 'otro'
        self.method = 'newton'
        self.L = 0
        self.libraries = False
        self.num_points = 100
        self.degree = None
        
        self.cid_press = self.plt_points.figure.canvas.mpl_connect('button_press_event', self.press_event)
        self.cid_move = self.plt_points.figure.canvas.mpl_connect('motion_notify_event', self.move_event)
        self.cid_release = self.plt_points.figure.canvas.mpl_connect('button_release_event', self.release_event)
        self.cid_erease = self.plt_points.figure.canvas.mpl_connect('button_press_event', self.erease_event)
        
    def update_data(self, knots, method, L, libraries, num_points, degree):
        self.knots = knots
        self.method = method
        self.L = L
        self.libraries = libraries
        self.num_points = num_points
        self.degree = degree if degree != 0 else None
    
    def print_data(self):
        print self.knots, ",", self.method, ",", self.L, ",", self.libraries, ",", self.num_points, ",", self.degree  
          
    def update_curve(self):
<<<<<<< HEAD
        if self.knots == 'otros': 
            knots = np.linespace(0,1,len(self.polygon))
        else:
            knots = self.knots
=======
        if self.knots == "otro": knots = np.linspace(0,1,len(self.polygon)) 
        else: knots = self.knots

>>>>>>> f2d9888bbe688006b78b7a12ec4ba19a10de1850
        p = polynomial_curve_fitting(np.array(self.polygon), knots, self.method, self.num_points, libraries=self.libraries, L=self.L, degree=self.degree)
        self.plt_curve.set_data(p[:,0],p[:,1])
        self.plt_points.set_data(self.xs, self.ys)
        self.plt_points.figure.canvas.draw()
        # self.print_data()
    
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
    root.title("Polynomial Curve Fitting")
    
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
    R2=Radiobutton(knots, text="Other    ", variable=vn, value=False)
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
    D = StringVar()
    D.set("0")
    Dvalue = Entry(method, textvariable=D, width=3, justify=CENTER, state='disabled')
    def check():
        if vm.get() == 0:
            Lvalue.configure(state='normal')
            Dvalue.configure(state='normal')
        else:
            Lvalue.configure(state='disabled')
            Dvalue.configure(state='disabled')
    R3=Radiobutton(method, text="Newton Polynomial   ", variable=vm, value=True, command=check)
    R4=Radiobutton(method, text="Least Square Fitting,    L =", variable=vm, value=False, command=check)
    R3.pack(anchor=W, side=LEFT)
    R4.pack(anchor=W, side=LEFT)
    Lvalue.pack(side=LEFT)
    mlabel2 = Label(method, text="  degree =")
    mlabel2.pack(side=LEFT)
    Dvalue.pack(side=LEFT)
    
    #Se crea la zona para seleccionar si se usan librerias o no
    library = Frame(root)
    library.pack(anchor=W) 
    llabel = Label(library, text="LIBRARIES:           ")
    llabel.pack(side=LEFT, padx=2, pady=2)
    vl = BooleanVar()
    vl.set(False)
    R5=Radiobutton(library, text="FALSE   ", variable=vl, value=False, command=check)
    R6=Radiobutton(library, text="TRUE", variable=vl, value=True, command=check)
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
        if vn.get() == 1:
            knots = 'chebyshev'
        else:
            knots = 'otro'
        #nk = int(Nvalue.get())
        if vm.get() == 0:
            method = 'least_squares'
            L = float(Lvalue.get())
            degree = int(Dvalue.get())
        else:
            method = 'newton'
            L = 0
            degree = None
        libraries = vl.get()
        num_points = int(Pvalue.get())
        linebuilder.update_data(knots, method, L, libraries, num_points, degree)
        linebuilder.update_curve()
    button = Button(buttons, text="UPDATE", command=update)
    button.pack()
    
    #Se añade un pequeño tutorial:
<<<<<<< HEAD
    T = Text(root, height=10, width=60)
    T.pack()
    T.insert(END, "INSTRUCTIONS: Para introducir el polinomio, se presiona con el botón izquierdo del ratón sobre el punto que se desee añadir.")
=======
    T = Text(root, height = 3)
    T.pack()
    T.insert(END, 
        '''INSTRUCTIONS: Left click to introduce points and find the interpolation curve. You can drag points to move them or right click to erase them. Configure above options to see different curves.'''
    )
>>>>>>> f2d9888bbe688006b78b7a12ec4ba19a10de1850
        
    root.mainloop()

if __name__ == "__main__": sys.exit(main())