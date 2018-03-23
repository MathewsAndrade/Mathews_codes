import numpy as np
import matplotlib.pyplot as plt
import warnings
from matplotlib import pyplot, widgets

def WGS84():
    '''
    This function returns the following parameters defining the 
    reference elipsoid WGS84:
    a = semimajor axis [m]
    f = flattening
    GM = geocentric gravitational constant of the Earth 
         (including the atmosphere) [m**3/s**2]
    omega = angular velocity [rad/s]
    
    output:
    a, f, GM, omega
    '''
    a = 6378137.0
    f = 1.0/298.257223563
    GM = 3986004.418*(10**8)
    omega = 7292115*(10**-11)
    
    return a, f, GM, omega
##################################################################################################################################

def gamma_somigliana(a, f, GM, omega, phi):
    '''
    This function calculates the normal gravity by using
    the Somigliana's formula.
    
    input:
    a: float containing the semimajor axis [m]
    f: float containing the flattening
    GM: float containing the geocentric gravitational constant 
        of the Earth (including the atmosphere) [m**3/s**2]
    omega: float containing the angular velocity [rad/s]
    phi: array containing the geodetic latitudes [degree]
    
    output:
    gamma: array containing the values of normal gravity
           on the surface of the elipsoid for each geodetic
           latitude [mGal]
    '''
    b = a*(1.0-f)
    a2 = a**2
    b2 = b**2
    E = np.sqrt(a2 - b2)
    elinha = E/b
    bE = b/E
    Eb = E/b
    atg = np.arctan(Eb)
    q0 = 0.5*((1+3*(bE**2))*atg - (3*bE))
    q0linha = 3.0*(1+(bE**2))*(1-(bE*atg)) - 1
    m = (omega**2)*(a2)*b/GM
    aux = elinha*q0linha/q0
    gammaa = (GM/(a*b))*(1-m-(m/6.0)*aux)
    gammab = (GM/a2)*(1+(m/3.0)*aux)
    aux = np.deg2rad(phi)
    s2 = np.sin(aux)**2
    c2 = np.cos(aux)**2
    # the 10**5 converts from m/s**2 to mGal
    gamma = (10**5)*((a*gammaa*c2) + (b*gammab*s2))/np.sqrt((a2*c2) + (b2*s2))
    return gamma
##################################################################################################################################
    
def gamma_closedform(a, f, GM, omega, phi, h):
    '''
    This function calculates the normal gravity by using
    a closed-form formula.
    
    input:
    a: float containing the semimajor axis [m]
    f: float containing the flattening
    GM: float containing the geocentric gravitational constant 
        of the Earth (including the atmosphere) [m**3/s**-2]
    omega: float containing the angular velocity [rad/s]
    phi: array containing the geodetic latitudes [degree]
    h: array containing the normal heights [m]
    
    output:
    gamma: array containing the values of normal gravity
           on the surface of the elipsoid for each geodetic
           latitude [mGal]
    '''
    b = a*(1.0-f)
    a2 = a**2
    b2 = b**2
    E = np.sqrt(a2 - b2)
    E2 = E**2
    bE = b/E
    Eb = E/b
    atanEb = np.arctan(Eb)
    phirad = np.deg2rad(phi)
    tanphi = np.tan(phirad)
    cosphi = np.cos(phirad)
    sinphi = np.sin(phirad)
    beta = np.arctan(b*tanphi/a)
    sinbeta = np.sin(beta)
    cosbeta = np.cos(beta)
    zl = b*sinbeta+h*sinphi
    rl = a*cosbeta+h*cosphi
    zl2 = zl**2
    rl2 = rl**2
    dll2 = rl2-zl2
    rll2 = rl2+zl2
    D = dll2/E2
    R = rll2/E2
    cosbetal = np.sqrt(0.5*(1+R) - np.sqrt(0.25*(1+R**2) - 0.5*D))
    cosbetal2 = cosbetal**2
    sinbetal2 = 1-cosbetal2
    bl = np.sqrt(rll2 - E2*cosbetal2)
    bl2 = bl**2
    blE = bl/E
    Ebl = E/bl
    atanEbl = np.arctan(Ebl)
    q0 = 0.5*((1+3*(bE**2))*atanEb - (3*bE))
    q0l = 3.0*(1+(blE**2))*(1-(blE*atanEbl)) - 1
    W = np.sqrt((bl2+E2*sinbetal2)/(bl2+E2))

    gamma = GM/(bl2+E2) - cosbetal2*bl*omega**2
    gamma += (((omega**2)*a2*E*q0l)/((bl2+E2)*q0))*(0.5*sinbetal2 - 1./6.)
    # the 10**5 converts from m/s**2 to mGal
    gamma = (10**5)*gamma/W
    
    return gamma
    
def estatistica(dado, unidade=None):
    '''
    Calcula os valores minimo, medio, maximo e 
    diferenca entre o maximo e o minimo de um
    dado.
    
    input
    
    dado: numpy array 1D - vetor de dados.
    unidade: string - unidade dos dados.
    
    output
    
    minimo: float - minimo dos dados.
    media: float - media dos dados.
    maximo: float - maximo dos dados.
    variacao: float - diferenca entre o maximo e o
              minimo dos dados.
    '''
    
    assert dado.size > 1, 'O vetor de dados deve ter mais de um elemento'
    
    minimo = np.min(dado)
    maximo = np.max(dado)
    media = np.mean(dado)
    variacao = maximo - minimo
    
    if (unidade != None):
        print ('     min.: %15.5f %s' % (minimo, unidade) )
        print ('    media: %15.5f %s' % (media, unidade) )
        print ('     max.: %15.5f %s' % (maximo, unidade) )
        print ('var. max.: %15.5f %s' % (variacao, unidade) )
    else:
        print ('     min.: %15.5f' % (minimo) )
        print ('    media: %15.5f' % (media) )
        print ('     max.: %15.5f' % (maximo) )
        print ('var. max.: %15.5f' % (variacao) )
        
    return minimo, media, maximo, variacao
##################################################################################################################################
    
def plota_mapa(projecao, x, y, dado, area, unidade, titulo, cores, tamanho,
               delta, perfis=None, estados=None, escala=None, eixos=None):
    '''
    Plota um mapa dos dados "dado", com coordenadas "x" e 
    "y" referidas a uma determinada projecao cartografica 
    "projecao". "unidade" e "titulo" sao, respectivamente,
    a unidade dos dados e o titulo do mapa.
    
    input
    
    projecao: objeto do mpl_toolkits.basemap.Basemap - 
              projecao cartografica.
    x: numpy array 1D - vetor com as coordenadas x
       da projecao.
    y: numpy array 1D - vetor com as coordenadas y
       da projecao.
    dado: numpy array 1D - vetor com os dados a serem
          plotados.
    area: list - lista com os valores minimo e maximo da
          longitude e minimo e maximo da latitude, em
          graus.
    unidade: string - unidade dos dados a serem plotados.
    titulo: string -  titulo do mapa.
    cores: codigo do colormaps_reference.py - esquema de 
           cores do mapa.
    tamanho: tuple - define os tamanhos tx e ty do mapa ao longo,
             repectivamente, dos eixos x e y. Os valores tx e ty
             sao passados entre parenteses e separados por virgula
             (tx, ty).
    delta: float - intervalo, em graus, entre os meridiando e paralelos.
    perfis: numpy array 2D - matriz com as coordenadas x e y
            dos pontos iniciais e finais dos perfis a serem 
            analisados. As coordenadas x e y estao na primeira
            e segunda colunas, respectivamente. As primeiras
            duas linhas contem os dois pontos que formam o 
            primeiro perfil, as proximas duas contem os pontos
            que formam o segundo perfil e assim sucessivamente.
    estados: boolean - se for igual a True, plota o contorno dos estados.
    escala: boolean - se for igual a True, plota a escala do mapa.
    eixos: boolean - se for igual a True, plota os eixos do mapa.
    
    output
    
    mapa: string - codigo de uma matplotlib.figure.Figure.
    '''

    dado_min = np.min(dado)
    dado_max = np.max(dado)
    
    #Esquema da escala de cores
    if (dado_min*dado_max < 0.):
        ranges = np.max(np.abs([dado_min, dado_max]))
        ranges_0 = 0.
    else:
        ranges = 0.5*(dado_max - dado_min)
        ranges_0 = 0.5*(dado_max + dado_min)
    
    longitude_central = 0.5*(area[1] + area[0])
    latitude_central = 0.5*(area[3] + area[2])
    
    x_max = np.max(x)*0.001 # valor maximo de x em km
    x_min = np.min(x)*0.001 # valor minimo de x em km
    
    if escala == True:
    
        #Valor em km a ser representado na escala
        #Este valor foi estabelecido como aproximadamente 
        #40 porcento da variacao maxima em x
        comprimento_escala = np.floor(0.4*(x_max - x_min)/100.)*100.
    
        #Posicao do centro da escala em coordenadas geodesicas
        longitude_escala = area[1] - 0.25*(area[1] - area[0])
        latitude_escala = area[2] + 0.05*(area[3] - area[2])
    
    x_min = np.min(x)
    x_max = np.max(x)
    y_min = np.min(y)
    y_max = np.max(y)
    
    plt.figure(figsize=tamanho)
    plt.title(titulo, fontsize=18, y=1.05)
    projecao.contourf(x, y, dado, 100, tri=True, cmap=plt.get_cmap(cores),
                      vmin = -ranges + ranges_0, vmax = ranges + ranges_0)
    plt.colorbar(orientation='horizontal', pad=0.04, aspect=50, 
                 shrink=0.7).set_label(unidade, fontsize=18)
    projecao.drawcoastlines()
    if (np.ceil(area[2]) == area[2]):
        parallels = np.arange(np.ceil(area[2]) + 1., area[3], delta)
    else:
        parallels = np.arange(np.ceil(area[2]), area[3], delta)
    if (np.ceil(area[0]) == area[0]):
        meridians = np.arange(np.ceil(area[0]) + 1., area[1], delta)
    else:
        meridians = np.arange(np.ceil(area[0]), area[1], delta)
    if eixos == True:
        projecao.drawparallels(parallels, labels=[1,1,0,0])
        projecao.drawmeridians(meridians, labels=[0,0,1,1])
    else:
        projecao.drawparallels(parallels)
        projecao.drawmeridians(meridians)
    if estados == True:
        projecao.drawstates()
    if perfis.all() != None: #temperfil == True:
        for i in range(0,perfis.shape[0],2):
            projecao.plot(perfis[i:i+2,0], perfis[i:i+2,1], 'o-k', linewidth=2)
    if escala == True:
        projecao.drawmapscale(longitude_escala, latitude_escala,
                            longitude_central, latitude_central,
                            length=comprimento_escala, barstyle='fancy')    
    plt.show()
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

