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



class FastIntersectionBezierCore(object):

    _EPSILON = 0.01
    
    @staticmethod
    def intersect_bbox_check(pA, pB, epsilon = _EPSILON):
        AX, AY = np.amax(pA,axis=0); Ax, Ay = np.amin(pA, axis=0)
        BX, BY = np.amax(pB,axis=0); Bx, By =  np.amin(pB, axis=0)
        return Ax <= BX and AX >= Bx and Ay <= BY and AY >= By

    @staticmethod
    def division(poly, t = 0.5):
        divA, divB = np.empty(poly.shape), np.empty(poly.shape)
        divA[0], divB[0] = poly[0], poly[-1]

        curr = np.copy(poly)
        for j in xrange(1,len(poly)): 
            np.add(curr,np.roll(curr,-1,axis=0), curr)
            curr = np.multiply(curr,t)[:-1]
            # np.multiply(curr,t,curr) # curr *= t

            divA[j], divB[j] = curr[0], curr[-1]

        return (divA, divB)

    @staticmethod
    def seg_intersect(a1,a2, b1,b2) :
        da = a2-a1; da[0], da[1] = da[1], -da[0]
        db = b2-b1
        return (np.dot(da,a1-b1)/np.dot(da,db))*db + b1
        # return db*np.dot(da,a1-b1)/np.dot(da,db) + b1



class IntersectionBezier(FastIntersectionBezierCore):
    _EPSILON = FastIntersectionBezierCore._EPSILON
    
    
    @staticmethod
    def _delta22(poly):
        a =  np.roll(poly,-1,0)
        np.subtract(poly,a,a)

        return np.subtract(a,np.roll(a,-1,0),a)[:-2]


    def __call__(self, polyA, polyB, epsilon = _EPSILON):
        '''
        args:
        - polyA: shape = (n,2)
        - polyB: shape = (m,2)
        
        kwargs:
        - epsilon: zero threashold

        return: k intersection points numpy array with shape = (k,2)
        '''
        self.lastPolyA, self.lastPolyB = polyA, polyB

        if self.intersect_bbox_check(polyA,polyB):
            deltaA = self._delta22(polyA)
            deltaB = self._delta22(polyB)

            DA, DB = np.amax(np.absolute(deltaA)), np.amax(np.absolute(deltaB))
            n, m, eps = len(polyA), len(polyB), epsilon*8

            if n*(n-1)*DA > eps: polyA, (divA,divB) = polyB, self.division(polyA)
            elif m*(m-1)*DB > eps: polyA, (divA,divB) = polyA, self.division(polyB)
            else: return [self.seg_intersect(polyA[0], polyA[-1],polyB[0], polyB[-1])] 

            ret = self(polyA,divA, epsilon = epsilon)
            ret.extend(self(polyA,divB, epsilon = epsilon))

            self.intersectionPoints = ret
            return ret
        self.intersectionPoints = []
        return []


    def _plot(self, polyA, k = 3):
        '''
        Plots the curve with k steps

        kwargs:
        - k: the number of subdivision steps
        - figure: the matplotlib figure to plot
        '''
        if len(polyA) == 0: return []
        if k == 0: 
            return polyA
        
        divA, divB = self.division(polyA)

        r = self._plot(divA, k - 1)
        return np.concatenate((r,self._plot(divB[::-1], k - 1)))




class DrawCurves(object):
    def __init__(self, pointsA, pointsB):
        self.plt_pointsA, self.plt_pointsB = pointsA, pointsB
        self.xsA, self.ysA = list(self.plt_pointsA.get_xdata()), list(self.plt_pointsA.get_ydata())
        self.xsB, self.ysB = list(self.plt_pointsB.get_xdata()), list(self.plt_pointsB.get_ydata())
        self.polygonA, self.polygonB = zip(self.xsA, self.ysA), zip(self.xsB, self.ysB)
        self.plt_curveA, = pointsA.axes.plot([], [], 'r')
        self.plt_curveB, = pointsB.axes.plot([], [], 'b')
        self.plt_inters, = pointsA.axes.plot([], [], 'go')
        self.selected_point = None
        self.selected_polygon = True
        
        self.cid_press = self.plt_pointsA.figure.canvas.mpl_connect('button_press_event', self.press_event)
        self.cid_move = self.plt_pointsA.figure.canvas.mpl_connect('motion_notify_event', self.move_event)
        self.cid_release = self.plt_pointsA.figure.canvas.mpl_connect('button_release_event', self.release_event)
        self.cid_erease = self.plt_pointsA.figure.canvas.mpl_connect('button_press_event', self.erease_event)

        self.intersector = IntersectionBezier()

    def update_Bezier(self):
        # añadir k
       pA = self.intersector._plot(np.array(self.polygonA), k = 3)
       self.plt_curveA.set_data([x for x,_ in pA], [y for _,y in pA])

       pB = self.intersector._plot(np.array(self.polygonB), k = 3)
       self.plt_curveB.set_data([x for x,_ in pB], [y for _,y in pB])

    
    def update_curve(self):
        "update the drawn curve"
        
        self.update_Bezier()
        self.plt_pointsA.set_data(self.xsA, self.ysA)
        self.plt_pointsA.figure.canvas.draw()
        self.plt_pointsB.set_data(self.xsB, self.ysB)
        self.plt_pointsB.figure.canvas.draw()


    def update_intersection(self):
        # añadir epsilon
        intr = self.intersector(
            np.array(self.polygonA), 
            np.array(self.polygonB),
            epsilon = 0.01)
        self.plt_inters.set_data([x for x,_ in intr], [y for _,y in intr])
        self.plt_inters.figure.canvas.draw()
    
    def new_point(self, event):
        if self.selected_polygon:
            self.xsA.append(event.xdata)
            self.ysA.append(event.ydata)
            self.polygonA.append([event.xdata,event.ydata])
        else:
            self.xsB.append(event.xdata)
            self.ysB.append(event.ydata)
            self.polygonB.append([event.xdata,event.ydata])
        print self.polygonA, self.polygonB

        self.update_curve()
        self.update_intersection()
         
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
        self.update_intersection() 
        
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
            self.update_intersection()
            

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
        linebuilder.update_curve()
        linebuilder.update_intersection()
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