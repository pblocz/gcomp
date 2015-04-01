# -*- coding: utf-8 -*-
from Tkinter import * 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.special import binom
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D 
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backend_bases import key_press_handler

class BezierPolicy(object):
    def __init__(self, polygon):
        self.polygon = np.array(polygon)
        self.N = len(self.polygon)
    def __call__(self, npoints=100, *args, **kwargs): 
        return self.compute(np.linspace(0, 1, npoints), 
                            *args, **kwargs)

class BernsteinPolicy(BezierPolicy):
    def _berstein(self,t):
        res = np.zeros((self.N+1, len(t),))
        for i in range(self.N):        
            res[i, :] = binom(self.N - 1 , i) * t**i * (1 - t)**(self.N - i -1)
        return res
    def compute(self, t):
        bernstein = self._berstein(t)
        curve_x = sum(self.polygon[i,0] * bernstein[i, :] \
                      for i in range(self.N))
        curve_y = sum(self.polygon[i,1] * bernstein[i, :] \
                      for i in range(self.N))
        return np.vstack((curve_x, curve_y,)).T

class DrawCurves(object):
    def __init__(self, pointsA, pointsB):
        self.plt_pointsA, self.plt_pointsB = pointsA, pointsB
        self.xsA, self.ysA = list(self.plt_pointsA.get_xdata()), list(self.plt_pointsA.get_ydata())
        self.xsB, self.ysB = list(self.plt_pointsB.get_xdata()), list(self.plt_pointsB.get_ydata())
        self.polygonA, self.polygonB = zip(self.xsA, self.ysA), zip(self.xsB, self.ysB)
        self.plt_curveA, = pointsA.axes.plot([], [], 'r')
        self.plt_curveB, = pointsB.axes.plot([], [], 'b')
        self.selected_point = None
        self.selected_polygon = True
        
        self.cid_press = self.plt_pointsA.figure.canvas.mpl_connect('button_press_event', self.press_event)
        self.cid_move = self.plt_pointsA.figure.canvas.mpl_connect('motion_notify_event', self.move_event)
        self.cid_release = self.plt_pointsA.figure.canvas.mpl_connect('button_release_event', self.release_event)
        self.cid_erease = self.plt_pointsA.figure.canvas.mpl_connect('button_press_event', self.erease_event)
    
    def update_Bezier(self):
        curve = BernsteinPolicy(self.polygonA)(1000)
        self.plt_curveA.set_data(curve[:,0], curve[:,1])
        curve = BernsteinPolicy(self.polygonB)(1000)
        self.plt_curveB.set_data(curve[:,0], curve[:,1])
    
    def update_curve(self):
        "update the drawn curve"
        self.update_Bezier()
        self.plt_pointsA.set_data(self.xsA, self.ysA)
        self.plt_pointsA.figure.canvas.draw()
        self.plt_pointsB.set_data(self.xsB, self.ysB)
        self.plt_pointsB.figure.canvas.draw()
    
    def new_point(self, event):
        if self.selected_polygon:
            self.xsA.append(event.xdata)
            self.ysA.append(event.ydata)
            self.polygonA.append([event.xdata,event.ydata])
        else:
            self.xsB.append(event.xdata)
            self.ysB.append(event.ydata)
            self.polygonB.append([event.xdata,event.ydata])
        self.update_curve()
         
    def press_event(self, event):
        if (event.inaxes != self.plt_pointsA.axes and event.inaxes != self.plt_pointsB.axes) or event.button != 1: 
            return
        if self.selected_polygon:
            state, indices = self.plt_pointsA.contains(event)
        else:
            state, indices = self.plt_pointsB.contains(event)                
        if not state:
            self.new_point(event)
        else:
            self.selected_point = indices['ind'][0]
            # print self.selected_point
            
    def move_event(self, event): 
        if (event.inaxes != self.plt_pointsA.axes and event.inaxes != self.plt_pointsB.axes) or self.selected_point is None:
            return
        if self.selected_polygon:
            self.xsA[self.selected_point] = event.xdata
            self.ysA[self.selected_point] = event.ydata
            self.polygonA[self.selected_point] = [event.xdata,event.ydata]
        else:
            self.xsB[self.selected_point] = event.xdata
            self.ysB[self.selected_point] = event.ydata
            self.polygonB[self.selected_point] = [event.xdata,event.ydata]
        self.update_curve()     
        
    def release_event(self, event):  
        if event.inaxes != self.plt_pointsA.axes and event.inaxes != self.plt_pointsB.axes: 
            return
        self.selected_point = None
        
    def erease_event(self, event):
        if (event.inaxes != self.plt_pointsA.axes and event.inaxes != self.plt_pointsB.axes) or event.button != 3: 
            return
        if self.selected_polygon:
            state, indices = self.plt_pointsA.contains(event)
        else:
            state, indices = self.plt_pointsB.contains(event)
        if state:
            index = indices['ind'][0]
            if self.selected_polygon:
                del self.xsA[index]
                del self.ysA[index]
                del self.polygonA[index]
            else:
                del self.xsB[index]
                del self.ysB[index]
                del self.polygonB[index]
            self.update_curve()
            

def main(args=None):
    args = args or sys.argv
    
    #Se crea la ventana
    root = Tk()
    root.title("Intersección de curvas de Bezier")
    
    #Se crea la zona para los botones de seleccion de poligonos
    label = Label(root, text="Elige un polígono")
    label.pack()
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
    
    #Se crea la zona de seleccion de datos
    label = Label(root, text="Elige los datos")
    label.pack()
    data = Frame(root)
    data.pack()
    labelK = Label(data, text="K")
    labelK.pack(side=LEFT, padx=2, pady=2)
    dataK = Entry(data)
    dataK.pack(side=LEFT, padx=2, pady=2)
    labelK = Label(data, text="              epsilon")
    labelK.pack(side=LEFT, padx=2, pady=2)
    dataE = Entry(data)
    dataE.pack(side=LEFT, padx=2, pady=2)
        
    #Se generan los poligonos y sus curvas de Bezier
    pointsA, = ax.plot([], [], 'ro')
    pointsB, = ax.plot([], [], 'bo')   
    linebuilder = DrawCurves(pointsA, pointsB)
       
    #Una vez generadas las curvas, se generar las funciones que tendran los botones
    def functionA():
        linebuilder.selected_polygon = True
    def functionB():
        linebuilder.selected_polygon = False
    def actualizar():
        print "Actualizado"  ################################Esto falta por hacer
        
    #Una vez que existen las funciones, se crean los botones y se les asignan las acciones 
    buttonA = Button(buttons, text="Polygon A", command=functionA, fg="red" )
    buttonA.pack(side=LEFT, padx=2, pady=2)
    buttonB = Button(buttons, text="Polygon B", command=functionB, fg="blue" )
    buttonB.pack(side=LEFT, padx=2, pady=2)
    buttonAct = Button(root, text="Actualizar", command=actualizar)
    buttonAct.pack()
    
    root.mainloop()

if __name__ == "__main__": sys.exit(main())