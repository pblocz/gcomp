import numpy as np
import matplotlib.pyplot as plt
from scipy.special import binom
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
import math

def keyfun(pivot):
    pix, piy = pivot
    def key( (px, py,) ):
        return (math.atan2(py - piy, px - pix), (pix-px)**2 + (piy-py)**2,)
    return key

def ccw((ax,ay), (bx,by), (cx,cy)):
    vlx, vly = bx-ax, by-ay
    vcx, vcy = cx-ax, cy-ay
    cross = vcy*vlx - vcx*vly

    return cross > 0

def convexhull(points):
    m = min(points, key = lambda (x,y,): (y,-x,))
    points = sorted(filter(lambda x: x!=m, points),key = keyfun(m))

    res = [points[-1],m , points[0]]
    for p in points[1:]:
        while not ccw(res[-2], res[-1], p): res.pop()
        res.append(p)
    return res

class CurvaDeBezier:
    def __init__(self, polygon):
        self.polygon = polygon
        self.N = polygon.shape[0] - 1      
        
        self.num_points = 100
        self.t = np.linspace(0, 1, self.num_points) 
        
        self.curve_x = None
        self.curve_y = None
    
        self._bernstein = np.zeros((self.N + 1, self.num_points))
        self._compute_berstein()
        self.compute_curve()
        
    def _compute_berstein(self):
        for i in range(self.N + 1):        
            self._bernstein[i, :] = binom(self.N, i) * self.t**i * (1 - self.t)**(self.N - i)
        
    def compute_curve(self):
        self.curve_x = sum(self.polygon[i, 0] * self._bernstein[i, :] for i in range(self.N + 1))
        self.curve_y = sum(self.polygon[i, 1] * self._bernstein[i, :] for i in range(self.N + 1))
        
        print sum(self.polygon[i,:] * self._bernstein[i,:] for i in range(self.N+1))

    def plot_bezier(self):
        self.curve = Line2D(self.curve_x, self.curve_y)
        ax.add_line(self.curve)
        fig.canvas.draw()
        
        
if __name__ == '__main__':
    
    points = [[0, 0.5], [1, 1], [0,3], [2, 0.5]]
    # convexhull(points)

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect=1)    
    ax.set_xlim(-1, 3)    
    ax.set_ylim(-1, 4)    
    
    # points = [[0, 0], [1, 1], [2, 0.5]]
    px = [x for x,_ in points]; py = [y for _,y in points]
    P = np.array(points)
    bezier = CurvaDeBezier(P)

    # plot points + polygon + convex hull
    plt.plot(px,py, "ro", figure=fig)
    plt.plot(px,py, "r", figure=fig)
    points = convexhull(points)
    arg = (lambda points: ([x for x,_ in points], [y for _,y in points],))(points+[points[0]])
    plt.plot(*(arg+("g",)), figure=fig)
    
    bezier.plot_bezier()

    plt.show()
    

    
