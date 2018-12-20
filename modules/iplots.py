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

# global variables for counting the number of clicks in both axes plots:
def model_masses(area1, area2, background=None):
    
    """ Function to plant 2D point masses by clicking with the mouse. 
    Inputs: 
    * area1 = [xmin, xmax, zmin, zmax] : list with 2D Cartesian coordinate ranges(x,z).
    * area2 = [rhomin, rhomax, rhomin, rhomax] : list with density-constrast and depth ranges.

           
    Output : x,z,rho = lists with the picked values from the mouse clicking. The size is related to the number of clicks
    OBS: PAY ATTENTION TO THE NUMBER OF CLICKS IN BOTH PLOT AREAS. OTHERWISE THE REMAINING CODE WILL NOT WORK! 
    """
 
    # ------------ auxiliar functios to perform the clicking--------------------------:
    def draw_guide1(px, py):
        if len(x) != 0 or len(y) !=0:
            tmpline1.set_data([x[-1], px], [y[-1], py])
    # --------------------------------------------------------------------------------:
    def draw_guide2(px,py):
        if len(rho) !=0 or len(z) !=0 :
            tmpline2.set_data([rho[-1], px], [z[-1], py])
    # --------------------------------------------------------------------------------:
    def move(event):
        if event.inaxes != ax1 or event.inaxes != ax2:
            return 'plot area wrongly set up. Please, check.'       
        if event.inaxes == ax1:
            if picking[0]:
                draw_guide1(event.xdata, event.ydata)
                ax1.figure.canvas.draw()
        elif event.inaxes == ax2:
            if picking[0]:
                draw_guide2(event.xdata, event.ydata)
                ax2.figure.canvas.draw()
    
    # --------------------------------------------------------------------------------:
    def click(event):
        if event.inaxes == ax1:
            # count for click instances:
            if event.button ==1:
                click1.append(1.0)
            #---- append list with picked values -----:
                x.append(event.xdata)
                y.append(event.ydata)
                plotx.append(event.xdata)
                ploty.append(event.ydata)
            
            #-------------- plot data -------------------: 
            line1.set_color('k')
            line1.set_marker('o')
            line1.set_linestyle('None')
            line1.set_data(plotx, ploty)
            # -------- display the number of clicks in the subtitle --------:
            ax1.set_title('Number of clicks ='+ str( len(click1) ),fontsize =12, color = 'black')
            
            ax1.figure.canvas.draw()
            
        elif event.inaxes == ax2:
            # count for click instances:
            if event.button ==1:
                click2.append(1.0)
                
             #---- append list with picked values:
            rho.append(event.xdata)
            z.append(event.ydata)
            plotrho.append(event.xdata)
            plotz.append(event.ydata)
            
            #-------------- plot data -------------------: 
            line2.set_color('b')
            line2.set_marker('*')
            line2.set_linestyle('None')
            line2.set_data(plotrho, plotz)
            # -------- display the number of clicks in the subtitle --------:
            ax2.set_title('Number of clicks =' + str( len(click2) ), fontsize =12, color = 'blue')
            ax2.figure.canvas.draw()
             
   # --------------------------------------------------------------------------------:
    def erase(event):
        if event.inaxes == ax1:
            if event.key == 'e' and picking[0]:
                # count for click instances: 
                click1.pop()     
             #---- remove list with "unpicked" values:
                x.pop()
                y.pop()
                plotx.pop()
                ploty.pop()
                    
             #-------------- plot data -------------------:
            line1.set_data(plotx, ploty)
            line1.set_linestyle('None')
            draw_guide1(event.xdata, event.ydata)
             # -------- display the number of clicks in the subtitle --------:
            ax1.set_title('Number of clicks =' + str( len(click1) ), fontsize =12, color = 'black')
            ax1.figure.canvas.draw()
        
        elif event.inaxes == ax2:
             if event.key == 'e' and picking[0]:
                 # count for click instances: 
                click2.pop() 
               # s2 = str(click2)
                rho.pop()
                z.pop()
                plotrho.pop()
                plotz.pop()
                #-------------- plot data -------------------:
                line2.set_data(plotrho, plotz)
                line2.set_linestyle('None')
                draw_guide2(event.xdata, event.ydata)
                # -------- display the number of clicks in the subtitle --------:
                ax2.set_title('Number of clicks =' + str( len(click2) ), fontsize =12, color = 'blue')
                ax2.figure.canvas.draw()
    #---------------------------------------------------------------------------------:
    
    fig1, (ax1, ax2) = pyplot.subplots(1, 2)
    
    fig1.suptitle('Click for (x, z) coordinates and (rho, rho) values.'
                  'Press keybord < e > to erase undesired values. Close figure when done', fontsize =16)
    ax1.set_xlim(area1[0], area1[1])
    ax1.set_ylim(area1[2], area1[3])
    ax1.grid()
    ax1.set_xlabel(' Horizontal coordinate x(m)')
    ax1.set_ylabel(' Depth (m)')
    ax1.invert_yaxis()
    
    # check for optional input:
    if background !=None:
        xv = background[0]
        zv = background[1]
        ax1.plot(xv,zv,'-k')
        
    ax2.set_xlim(area2[0], area2[1])
    ax2.set_ylim(area2[2], area2[3])
    ax2.grid()

# --------- change the position of the tick label--------------:
    ax2.yaxis.tick_right()

# --------- invert axis for depth positive down --------------:
    ax2.set_xlabel('density')
    ax2.set_ylabel('density')
     
# --------- cursor use for better visualization ------------- :
    cursor1 = widgets.Cursor(ax1, useblit=True, color='black', linewidth=2)
    cursor2 = widgets.Cursor(ax2, useblit=True, color='blue', linewidth=2)

    # ----- for the case of new clicks -------:
    x = []
    y = []
    z = []
    
    plotx = []
    ploty = []

    rho = []
    plotz = []
    plotrho = []
    click1 = []
    click2 = []
   # ----------------- cleaning line object for plotting ------------------:
    line1, = ax1.plot([],[])
    tmpline1, = ax1.plot([],[])
    line2, = ax2.plot([],[])
    tmpline2, = ax2.plot([],[])
    
# ------------ Hack because Python 2 doesn't like nonlocal variables that change value -------:
# Lists it doesn't mind.
    picking = [True]
    fig1.canvas.mpl_connect('button_press_event', click )
    fig1.canvas.mpl_connect('key_press_event', erase )
    fig1.canvas.mpl_connect('motion_notify_event', move )
    pyplot.show()
    return x,y,rho





