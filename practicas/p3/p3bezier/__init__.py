import sys
from .policies import FastBernsteinPolicy

class Bezier:
    def __init__(self, polygon, policy = FastBernsteinPolicy):
        self.polygon = polygon        
        self.num_points = 100
        
        self.policy = policy(polygon)
        self.curve = self.policy(npoints = self.num_points)
    
        self.curve_x = self.curve[:,0]
        self.curve_y = self.curve[:,1]


    def plot_bezier(self):
        self.curve = Line2D(self.curve_x, self.curve_y)
        ax.add_line(self.curve)
        fig.canvas.draw()
        

def main(args=None):
    args = args or sys.argv

    return 0

if __name__ == "__main__": sys.exit(main())
