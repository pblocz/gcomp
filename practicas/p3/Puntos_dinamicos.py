import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from mis_bonitas_geodesicas import Geodesica

class DrawPoints(object):
    def __init__(self, fig, ax, geodesica):
        self.fig = fig
        self.ax = ax
        self.exists_touched_circle = False
        self.cid_press = fig.canvas.mpl_connect('button_press_event',self.on_press)
        self.cid_move = fig.canvas.mpl_connect('motion_notify_event',self.on_move)
        self.cid_release = fig.canvas.mpl_connect('button_release_event',self.on_release)
    
    def on_press(self, event):   
        for circle in self.ax.patches:
            contains, attr = circle.contains(event)
            if contains:
                self.touched_circle = circle
                self.exists_touched_circle = True
                self.pressed_event = event
                self.touched_x0, self.touched_y0 = circle.center
                return
        c = Circle((event.xdata,event.ydata), 0.5)
        self.ax.add_patch(c)
        self.fig.canvas.draw()
    
    def on_move(self, event):
        if self.exists_touched_circle:
            dx = event.xdata - self.pressed_event.xdata
            dy = event.ydata - self.pressed_event.ydata
            x0,y0 = self.touched_circle.center
            self.touched_circle.center = self.touched_x0 +dx, self.touched_y0 + dy
            #calculos (actualizan la curva)
            self.fig.canvas.draw()
            
    def on_release(self, event):
        self.exists_touched_circle = False
        return
            
if __name__ == '__main__':
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect=1)
    ax.set_xlim(-20, 20)
    ax.set_ylim(-20, 20)
    draw_points = DrawPoints(fig,ax)
    plt.show()
