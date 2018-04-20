import numpy as np
##################################################################################################################################
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

##################################################################################################################################

def grav2D_anom(xv,zv,rho,beta):
# COMPUTES THE 2D GRAVITY ANOMALY IN mGals OF A TWO DIMENSIONAL(2D) BODY OF ARBITRARY SHAPE WITH HYPERBOLIC
#! DENBITY CONTRAST
#! This code is based on Rao et al 1994 implementation
#! xv,zv = vertices of the arbitrary source (kilometers) 
#! rho = density contrast in g/cm3
#! beta = coeficient to control the hyperbolic decay with depth 
#! n = number of vertices of the arbitrary source (in clockwise direction)
# output: grav (1D array with gravity anomaly produced by the polygon)
    
    # Check inconstancy of array dimensions:
    # Stablishing some conditions
    if xv.shape != zv.shape:
        raise ValueError("All inputs must have same shape!")
    # get the number of elements of xv (number of vertices of the polygon)
    n = np.size(xv)
        
    # create new working arrays for the vertices of a closed polygon:
    x = np.zeros( (n+1,) )  
    z = np.zeros( (n+1,) ) 
    x[0:n] = xv
    z[0:n] = zv
    # Closed body: 
    x[n:n+1] = xv[0]
    z[n:n+1] = zv[0]
    gval = 0.0    #! gravity value to be added
    #! Loop over the vertices of the 2D polygon:
    for k in range(n):
        rk = np.sqrt(x[k]**2 + z[k]**2)
        if rk == 0.0: 
            raise ValueError("Identical vertices. Check coordinates")
            break
        rk1 = np.sqrt( x[k+1]**2 + z[k+1]**2 )
        if rk1 == 0.0: 
            raise ValueError("Identical vertices. Check coordinates")
            break
        #! vertices differences:
        dx  = x[k+1] - x[k]
        dz  = z[k+1] - z[k]
        den = np.sqrt(dx*dx + dz*dz)
        if den == 0.0: 
            raise ValueError("Identical vertices. Check coordinates")
            break
  
    # compute the hyperbolic density function (Litinsky, 1994):
        c  = dx/den
        s  = dz/den
        p1 = x[k] * s - z[k] * c
        p2 = beta**2 - 2.0 * beta * p1 * c + p1**2
        q1 = beta + z[k]
        q2 = beta + z[k+1]
        if z[k] == 0.0:
            ph1 = 1.570796 * ( x[k] / np.absolute( x[k] ) )
        else:
            ph1 = np.arctan2(x[k], z[k] ) 

        if z[k+1] == 0.0:
            ph2 = 1.570796 * ( x[k+1] / np.absolute( x[k+1] ) )
        else:
            ph2 = np.arctan2( x[k+1], z[k+1] )
        
       #! Compute the gravity effect:
        d1 = rk1 * q1
        d2 = rk  * q2
        a1 = np.log(d1 / d2)
        t1 = (a1 * s) / p2
        t2 = (beta - p1 * c) / (p1 * p2)
        t3 = ph1 - ph2
        t4 = t2 * t3
        dg = p1 * (t1 - t4) 
   
        gval += dg
    
    cte = 13.33333 * rho * beta**2
    #!grav2D_anom = gval * cte     ! gravity anomaly in mGal
    grav = gval * cte
    return grav
##################################################################################################################################
