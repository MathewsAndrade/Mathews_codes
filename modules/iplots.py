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
################################################################################################################################
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
    * ax : matplotlib Axes
        The figure to use for drawing the polygon.
        To get an Axes instace, just do::
            from matplotlib import pyplot
            ax = pyplot.figure().add_subplot(1,1,1)
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
    #fig = pyplot.figure()
    #ax = fig.add_subplot(1,1,1)
    axes.set_title('Click to pick points (ALWAYS CLOCKWISE!!!). Close fig when done.')
    if xy2ne:
        axes.set_xlim(area[2], area[3])
        axes.set_ylim(area[0], area[1])
    else:
        axes.set_xlim(area[0], area[1])
        axes.set_ylim(area[2], area[3])
    
    axes.grid()
    line, = axes.plot([],[])
    tmpline, = axes.plot([], [])
    line.figure.canvas.draw()
    x = []
    y = []
    plotx = []
    ploty = []
    pyplot.gca().invert_yaxis()
    axes.figure.canvas.draw()
    # Hack because Python 2 doesn't like nonlocal variables that change value.
    # Lists it doesn't mind.
    picking = [True]
    def draw_guide(px, py):
        if len(x) != 0:
            tmpline.set_data([x[-1], px], [y[-1], py])

    def move(event):
        if event.inaxes != ax:
            return 'plot area wrongly set up. Please, check.'
        if picking[0]:
            draw_guide(event.xdata, event.ydata)
            axes.figure.canvas.draw()

    def pick(event):
        if event.inaxes != axes:
            return 'plot area wrongly set up. Please, check.'
        x.append(event.xdata)
        y.append(event.ydata)
        plotx.append(event.xdata)
        ploty.append(event.ydata)
        line.set_color('r')
        line.set_marker('o')
        line.set_linestyle('-')
        line.set_data(plotx,ploty)
        axes.figure.canvas.draw()

    def move(event):
        if event.inaxes != axes:
            return 'plot area wrongly set up. Please, check.'
        if picking[0]:
            draw_guide(event.xdata, event.ydata)
            axes.figure.canvas.draw()
        
    def erase(event):
        if event.key == 'e' and picking[0]:
            x.pop()
            y.pop()
            plotx.pop()
            ploty.pop()
            line.set_data(plotx, ploty)
            draw_guide(event.xdata, event.ydata)
            axes.figure.canvas.draw()

    line.figure.canvas.mpl_connect('button_press_event', pick)
    line.figure.canvas.mpl_connect('key_press_event', erase)
    line.figure.canvas.mpl_connect('motion_notify_event', move)
    pyplot.show()
    #if xy2ne:
    #    points = np.transpose([y, x])
    #else:
    #    points = np.transpose([x, y])
    return x,y
##################################################################################################################################
