import warnings
import numpy
from matplotlib import pyplot, widgets
# Quick hack so that the docs can build using the mocks for readthedocs
# Ideal would be to log an error message saying that functions from pyplot
# were not imported
try:
    from matplotlib.pyplot import *
except:
    pass
##################################################################################################################################
def draw_polygon(area, axes, style='-', marker='o', color='k', width=2,
                 alpha=0.5, xy2ne=False):
    """
    Draw a polygon by clicking with the mouse.
    INSTRUCTIONS:
    * Left click to pick the edges of the polygon;
    * Draw edges CLOCKWISE;
    * Press 'e' to erase the last edge;
    * Right click to close the polygon;
    * Close the figure window to finish;
    Parameters:
    * area : list = [x1, x2, y1, y2]
        Borders of the area containing the polygon
    * axes : matplotlib Axes
        The figure to use for drawing the polygon.
        To get an Axes instace, just do::
            from matplotlib import pyplot
            axes = pyplot.figure().add_subplot(1,1,1)
        You can plot things to ``axes`` before calling this function so that
        they'll appear on the background.
    * style : str
        Line style (as in matplotlib.pyplot.plot)
    * marker : str
        Style of the point markers (as in matplotlib.pyplot.plot)
    * color : str
        Line color (as in matplotlib.pyplot.plot)
    * width : float
        The line width (as in matplotlib.pyplot.plot)
    * alpha : float
        Transparency of the fill of the polygon. 0 for transparent, 1 for
        opaque (fills the polygon once done drawing)
    * xy2ne : True or False
        If True, will exchange the x and y axis so that x points north.
        Use this when drawing on a map viewed from above. If the y-axis of the
        plot is supposed to be z (depth), then use ``xy2ne=False``.
    Returns:
    * edges : list of lists
        List of ``[x, y]`` pairs with the edges of the polygon
    """
    axes.set_title("Click to draw polygon. Right click when done.")
    if xy2ne:
        axes.set_xlim(area[2], area[3])
        axes.set_ylim(area[0], area[1])
    else:
        axes.set_xlim(area[0], area[1])
        axes.set_ylim(area[2], area[3])
    # start with an empty line
    line, = axes.plot([], [], marker=marker, linestyle=style, color=color,
                      linewidth=width)
    tmpline, = axes.plot([], [], marker=marker, linestyle=style, color=color,
                         linewidth=width)
    draw = axes.figure.canvas.draw
    x = []
    y = []
    plotx = []
    ploty = []
    # Hack because Python 2 doesn't like nonlocal variables that change value.
    # Lists it doesn't mind.
    picking = [True]

    def draw_guide(px, py):
        if len(x) != 0:
            tmpline.set_data([x[-1], px], [y[-1], py])

    def move(event):
        if event.inaxes != axes:
            return 0
        if picking[0]:
            draw_guide(event.xdata, event.ydata)
            draw()

    def pick(event):
        if event.inaxes != axes:
            return 0
        if event.button == 1 and picking[0]:
            x.append(event.xdata)
            y.append(event.ydata)
            plotx.append(event.xdata)
            ploty.append(event.ydata)
        if event.button == 3 or event.button == 2 and picking[0]:
            if len(x) < 3:
                axes.set_title("Need at least 3 points to make a polygon")
            else:
                picking[0] = False
                axes.set_title("Done! You can close the window now.")
                plotx.append(x[0])
                ploty.append(y[0])
                tmpline.set_data([], [])
                axes.fill(plotx, ploty, color=color, alpha=alpha)
        line.set_data(plotx, ploty)
        draw()

    def erase(event):
        if event.key == 'e' and picking[0]:
            x.pop()
            y.pop()
            plotx.pop()
            ploty.pop()
            line.set_data(plotx, ploty)
            draw_guide(event.xdata, event.ydata)
            draw()
    line.figure.canvas.mpl_connect('button_press_event', pick)
    line.figure.canvas.mpl_connect('key_press_event', erase)
    line.figure.canvas.mpl_connect('motion_notify_event', move)
    pyplot.show()
    if len(x) < 3:
        raise ValueError("Need at least 3 points to make a polygon")
    if xy2ne:
        verts = numpy.transpose([y, x])
    else:
        verts = numpy.transpose([x, y])
    return verts
##################################################################################################################################
def pick_points(area, axes, marker='o', color='k', size=8, xy2ne=False):
    """
    Get the coordinates of points by clicking with the mouse.
    INSTRUCTIONS:
    * Left click to pick the points;
    * Press 'e' to erase the last point picked;
    * Close the figure window to finish;
    Parameters:
    * area : list = [x1, x2, y1, y2]
        Borders of the area containing the points
    * axes : matplotlib Axes
        The figure to use for drawing the polygon.
        To get an Axes instace, just do::
            from matplotlib import pyplot
            axes = pyplot.figure().add_subplot(1,1,1)
        You can plot things to ``axes`` before calling this function so that
        they'll appear on the background.
    * marker : str
        Style of the point markers (as in matplotlib.pyplot.plot)
    * color : str
        Line color (as in matplotlib.pyplot.plot)
    * size : float
        Marker size (as in matplotlib.pyplot.plot)
    * xy2ne : True or False
        If True, will exchange the x and y axis so that x points north.
        Use this when drawing on a map viewed from above. If the y-axis of the
        plot is supposed to be z (depth), then use ``xy2ne=False``.
    Returns:
    * points : list of lists
        List of ``[x, y]`` coordinates of the points
    """
   axes.set_title("Click to pick points. Close window when done.")
   if xy2ne:
      axes.set_xlim(area[2], area[3])
      axes.set_ylim(area[0], area[1])
   else:
      axes.set_xlim(area[0], area[1])
      axes.set_ylim(area[2], area[3])
    # start with an empty set
   line, = axes.plot([], [], marker=marker, color=color, markersize=size)
   line.figure.canvas.draw()
   x = []
   y = []
   plotx = []
   ploty = []
   # Hack because Python 2 doesn't like nonlocal variables that change value.
   # Lists it doesn't mind.
   picking = [True]
   def pick(event):
      if event.inaxes != axes:
         return 0
      if event.button == 1 and picking[0]:
         x.append(event.xdata)
         y.append(event.ydata)
         plotx.append(event.xdata)
         ploty.append(event.ydata)
         line.set_color(color)
         line.set_marker(marker)
         line.set_markersize(size)
         line.set_linestyle('')
         line.set_data(plotx, ploty)
   line.figure.canvas.draw()

   def erase(event):
      if event.key == 'e' and picking[0]:
         x.pop()
         y.pop()
         plotx.pop()
         ploty.pop()
         line.set_data(plotx, ploty)
   line.figure.canvas.draw()
   line.figure.canvas.mpl_connect('button_press_event', pick)
   line.figure.canvas.mpl_connect('key_press_event', erase)
   pyplot.show()
   if xy2ne:
      points = np.transpose([y, x])
   else:
      points = np.transpose([x, y])
   return points
##################################################################################################################################
