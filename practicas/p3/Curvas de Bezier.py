import numpy as np
import matplotlib.pyplot as plt
from scipy.special import binom
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D


class CurvaDeBezier(object):
    def __init__(self, polygon):
        self.polygon = polygon
        self.N = polygon.shape[0]
        self.num_points = 100
        self.t = np.linspace(0,1,self.num_points)
        self.curve_x = None
        self.curve_y = None
        self._berstein = np.zeros((self.N,self.num_points))
        self._compute_berstein()
        self._compute_curve()
        
    def _compute_berstein(self):
        for i in range(self.N):
            self._berstein[i,:] = binom(self.N,i)*self.t**i * (1-self.t)**(self.N-i)
    
    def _compute_curve(self):
        self.curve_x = sum(self.polygon[i,0] * self._berstein[i,:] for i in range(self.N))
        self.curve_y = sum(self.polygon[i,1] * self._berstein[i,:] for i in range(self.N))
    
    def plot_bezier(self):
        self.curve = Line2D(self.curve_x, self.curve_y)
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, aspect=1)
        
        self.ax.add_line(self.curve)
        self.fig.canvas.draw()     
        self.fig.show()
        
p = np.array([[0,0],[1,1],[2,0.5]])
CurvaDeBezier(p).plot_bezier()
