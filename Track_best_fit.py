# Python library for computing the best fit center of a satellite track. 
# 
# There are 5 possible functions that can be used:
#
# 1-With the function "track_center_superellipse(x, y)" the best fit is made with a superellipse with n=6 (slow method) 
#
# 2-With the function "track_center(x, y)" the median of the 4 independent sides of the track is first calculated to smooth the perimeter
# and then find the center of the trace thus obtained.
#
# 3-With the function "track_center2(x, y)" the outliers of the upper and lower edges of the track placed in the canonical system 
# are first eliminated. Then, as in method 2, the median of the 4 independent sides of the new track is calculated to smooth the perimeter
# finally the center of the trace thus obtained is found.
#
# 4-With the function "track_center_median(x, y)", the median of the vectors x and y of the raw perimeter of the track is calculated.
#
# 5-With the function "track_center_mean(x, y)", the mean of the vectors x and y of the raw perimeter of the track is calculated
#
# With satellite tracks the best results are obtained with 2 and 3.
#
# The other functions of this library support those listed above.
#
# Input:
# Vectors of the x and y coordinates in pixels of the track to fit
#
# Output:
# Best track center in pixels
#
# Albino Carbognani, INAF-OAS
# Versione del 4 agosto 2022

import numpy as np
import matplotlib.pyplot as plt
import statistics
from scipy.signal import savgol_filter

#==================================================================================

def fit_ellipse(x, y):

    """
    Fit the coefficients a,b,c,d,e,f, representing an ellipse described by
    the formula F(x,y) = ax^2 + bxy + cy^2 + dx + ey + f = 0 to the provided
    arrays of data points x=[x1, x2, ..., xn] and y=[y1, y2, ..., yn].

    Based on the algorithm of Halir and Flusser, "Numerically stable direct
    least squares fitting of ellipses'.
    https://scipython.com/blog/direct-linear-least-squares-fitting-of-an-ellipse/
    """

    D1 = np.vstack([x**2, x*y, y**2]).T
    D2 = np.vstack([x, y, np.ones(len(x))]).T
    S1 = D1.T @ D1
    S2 = D1.T @ D2
    S3 = D2.T @ D2
    T = -np.linalg.inv(S3) @ S2.T
    M = S1 + S2 @ T
    C = np.array(((0, 0, 2), (0, -1, 0), (2, 0, 0)), dtype=float)
    M = np.linalg.inv(C) @ M
    eigval, eigvec = np.linalg.eig(M)
    con = 4 * eigvec[0]* eigvec[2] - eigvec[1]**2
    ak = eigvec[:, np.nonzero(con > 0)[0]]
    return np.concatenate((ak, T @ ak)).ravel()

#==================================================================================

def cart_to_pol(coeffs):

    """
    Convert the cartesian conic coefficients, (a, b, c, d, e, f), to the
    ellipse parameters, where F(x, y) = ax^2 + bxy + cy^2 + dx + ey + f = 0.
    The returned parameters are x0, y0, ap, bp, e, phi, where (x0, y0) is the
    ellipse centre; (ap, bp) are the semi-major and semi-minor axes,
    respectively; e is the eccentricity; and phi is the rotation of the semi-
    major axis from the x-axis.
    https://scipython.com/blog/direct-linear-least-squares-fitting-of-an-ellipse/
    """

    # We use the formulas from https://mathworld.wolfram.com/Ellipse.html
    # which assumes a cartesian form ax^2 + 2bxy + cy^2 + 2dx + 2fy + g = 0.
    # Therefore, rename and scale b, d and f appropriately.
    a = coeffs[0]
    b = coeffs[1] / 2
    c = coeffs[2]
    d = coeffs[3] / 2
    f = coeffs[4] / 2
    g = coeffs[5]

    den = b**2 - a*c
    if den > 0:
        raise ValueError('coeffs do not represent an ellipse: b^2 - 4ac must'
                         ' be negative!')

    # The location of the ellipse centre.
    x0, y0 = (c*d - b*f) / den, (a*f - b*d) / den

    num = 2 * (a*f**2 + c*d**2 + g*b**2 - 2*b*d*f - a*c*g)
    fac = np.sqrt((a - c)**2 + 4*b**2)
    # The semi-major and semi-minor axis lengths (these are not sorted).
    ap = np.sqrt(num / den / (fac - a - c))
    bp = np.sqrt(num / den / (-fac - a - c))

    # Sort the semi-major and semi-minor axis lengths but keep track of
    # the original relative magnitudes of width and height.
    width_gt_height = True
    if ap < bp:
        width_gt_height = False
        ap, bp = bp, ap

    # The eccentricity.
    r = (bp/ap)**2
    if r > 1:
        r = 1/r
    e = np.sqrt(1 - r)

    # The angle of anticlockwise rotation of the major-axis from x-axis.
    if b == 0:
        phi = 0 if a < c else np.pi/2
    else:
        phi = np.arctan((2.*b) / (a - c)) / 2
        if a > c:
            phi += np.pi/2
    if not width_gt_height:
        # Ensure that phi is the angle to rotate to the semi-major axis.
        phi += np.pi/2
    phi = phi % np.pi

    return x0, y0, ap, bp, e, phi

#==================================================================================

def get_ellipse_pts(params, npts=200, tmin=0, tmax=2*np.pi):

    """
    Return npts points on the ellipse described by the params = x0, y0, ap,
    bp, e, phi for values of the parametric variable t between tmin and tmax.
    https://scipython.com/blog/direct-linear-least-squares-fitting-of-an-ellipse/
    """

    x0, y0, ap, bp, e, phi = params
    # A grid of the parametric variable, t.
    t = np.linspace(tmin, tmax, npts)
    x = x0 + ap * np.cos(t) * np.cos(phi) - bp * np.sin(t) * np.sin(phi)
    y = y0 + ap * np.cos(t) * np.sin(phi) + bp * np.sin(t) * np.cos(phi)

    # Ellipse's points
    return x, y

#==================================================================================

def get_superellipse_pts(params, npts=200, tmin=0, tmax=2*np.pi):

    """
    Return npts points on the superellipse described by the params = x0, y0, ap,
    bp, phi for values of the parametric variable t between tmin and tmax.
    A. Carbognani, INAF-OAS, July 27, 2022
    """

    x0, y0, ap, bp, phi = params

    # A grid of the parametric variable, t.
    t = np.linspace(tmin, tmax, npts)
    
    # Superellipse in canonical form with n=6
    x1=ap * np.sign(np.cos(t)) * (abs(np.cos(t)))**(1./3) 
    y1=bp * np.sign(np.sin(t)) * (abs(np.sin(t)))**(1./3)

    # Rotation and translation to get non-canonical superellipse
    x = x0 + x1 * np.cos(phi) - y1 * np.sin(phi)
    y = y0 + x1 * np.sin(phi) + y1 * np.cos(phi)  

    # Superellipse's points
    return x, y

#==================================================================================

def smooth_superellipse_pts(params, X, Y):

    """
    Smooth of the points X, Y using the superellipse n=6 described by the params = x0, y0, ap,
    bp, phi.
    A. Carbognani, INAF-OAS, July 27, 2022
    """

    ap, bp, phi = params

    # Put the best fit center of the superellipse 
    # to make them less sensitive to edge irregularities.
    x0=statistics.mean(X)
    y0=statistics.mean(Y)
    
    # Normal superellipse reduction
    X1=((X-x0)*np.cos(phi)+(Y-y0)*np.sin(phi))/ap;
    Y1=((Y-y0)*np.cos(phi)-(X-x0)*np.sin(phi))/bp;
    Q=(X1*X1*X1*X1*X1*X1 + Y1*Y1*Y1*Y1*Y1*Y1 - 1);
    j=0
    X2=[]
    Y2=[]
    for i in range(0, len(Q)):
        if Q[i] < 4.0:
           X2.insert(j, X[i])
           Y2.insert(j, Y[i])
           j=j+1 

    return X2, Y2

#==================================================================================

def smooth_borders_median(params, X, Y):

    """
    Compute track's best center using border's median algorithm.
    Input:
    x0, y0 = track's center, phi inclination angle (rad); 

    Output:
    X, Y = track's border
    A. Carbognani, INAF-OAS, July 29, 2022
    """

    x0, y0, phi = params

    # Reduction to canonical ellipse system (a=b=1)
    X1=((X-x0)*np.cos(phi)+(Y-y0)*np.sin(phi));
    Y1=((Y-y0)*np.cos(phi)-(X-x0)*np.sin(phi));
        
    # Separate track's borders
    Yup=[]
    Ydo=[]
    Xle=[]
    Xri=[]
    j=k=l=m=0
    MX1=0.96*max(X1)
    mX1=0.96*min(X1)

    for i in range(0, len(X1)):
        if X1[i] > mX1 and X1[i] < MX1 and Y1[i] > 0:
           Yup.insert(j, Y1[i])
           j=j+1
        if X1[i] > mX1 and X1[i] < MX1 and Y1[i] < 0:
           Ydo.insert(k, Y1[i])
           k=k+1
        if X1[i] >= min(X1) and X1[i] <= mX1:
           Xle.insert(l, X1[i])
           l=l+1
        if X1[i] >= MX1 and X1[i] <= max(X1):
           Xri.insert(l, X1[i])
           m=m+1

    # Compute borders's median values
    med_up=statistics.median(Yup)
    med_do=statistics.median(Ydo)
    med_le=statistics.median(Xle)
    med_ri=statistics.median(Xri)

    # Computer track's center in canonical system
    Xc=(med_le+med_ri)/2.0
    Yc=(med_up+med_do)/2.0
        
    # Compute median borders
    for i in range(0, len(X1)):
        if X1[i] > mX1 and X1[i] < MX1 and Y1[i] > 0:
           Y1[i]=med_up
        if X1[i] > mX1 and X1[i] < MX1 and Y1[i] < 0:
           Y1[i]=med_do
        if X1[i] >= min(X1) and X1[i] <= mX1:
           X1[i]=med_le
        if X1[i] >= MX1 and X1[i] <= max(X1):
           X1[i]=med_ri

    # Smooth borders rotation and translation to put in non-canonical system
    X2 = x0 + X1 * np.cos(phi) - Y1 * np.sin(phi)
    Y2 = y0 + X1 * np.sin(phi) + Y1 * np.cos(phi) 

    # Center rotation and translation to put in non-canonical system
    X0 = x0 + Xc * np.cos(phi) - Yc * np.sin(phi)
    Y0 = y0 + Xc * np.sin(phi) + Yc * np.cos(phi)

    return X2, Y2, X0, Y0

#==================================================================================

def smooth_borders_median2(params, X, Y):

    """
    Delete outliers in up and down track's border (in canonical system)
    Input:
    x0, y0 = track's center, phi inclination angle (rad);
 
    Output:
    Xn, Yn = No-outliers track's border
    A. Carbognani, INAF-OAS, Aug 03, 2022
    """

    x0, y0, phi = params

    # Reduction to canonical ellipse system (a=b=1)
    X1=((X-x0)*np.cos(phi)+(Y-y0)*np.sin(phi));
    Y1=((Y-y0)*np.cos(phi)-(X-x0)*np.sin(phi));
        
    # Separate track's borders
    Yup=[]
    Ydo=[]
    Xle=[]
    Xri=[]
    j=k=l=m=0
    MX1=0.96*max(X1)
    mX1=0.96*min(X1)

    for i in range(0, len(X1)):
        if X1[i] > mX1 and X1[i] < MX1 and Y1[i] > 0:
           Yup.insert(j, Y1[i])
           j=j+1
        if X1[i] > mX1 and X1[i] < MX1 and Y1[i] < 0:
           Ydo.insert(k, Y1[i])
           k=k+1
        if X1[i] >= min(X1) and X1[i] <= mX1:
           Xle.insert(l, X1[i])
           l=l+1
        if X1[i] >= MX1 and X1[i] <= max(X1):
           Xri.insert(l, X1[i])
           m=m+1

    # Compute borders's median values
    med_up=statistics.median(Yup)
    med_do=statistics.median(Ydo)
    med_le=statistics.median(Xle)
    med_ri=statistics.median(Xri)

    # Compute sigma for track's border up and down (in canonical system)
    sigma_up=[]
    sigma_do=[]
    j=k=0
    for i in range(0, len(X1)):
        if X1[i] > mX1 and X1[i] < MX1 and Y1[i] > 0:
            sigma_up.insert(j, (Y1[i]-med_up)*(Y1[i]-med_up))
            j=j+1
        if X1[i] > mX1 and X1[i] < MX1 and Y1[i] < 0:
            sigma_do.insert(k, (Y1[i]-med_do)*(Y1[i]-med_do))
            k=k+1

    Sigma_up=np.sqrt(sum(sigma_up)/len(sigma_up))
    Sigma_do=np.sqrt(sum(sigma_do)/len(sigma_do))
    
    # Erase outliers and save new track's border
    Yn=[]
    Xn=[]
    j=0

    for i in range(0, len(X1)):   
        if X1[i] >= min(X1) and X1[i] <= mX1:
           Xn.insert(j, X1[i])
           Yn.insert(j, Y1[i])
           j=j+1    

    for i in range(0, len(X1)):
        if X1[i] > mX1 and X1[i] < MX1 and Y1[i] > 0 and abs(Y1[i]-med_up) < 1.0*Sigma_up:
           Yn.insert(j, Y1[i])
           Xn.insert(j, X1[i])
           j=j+1

    for i in range(0, len(X1)):   
        if X1[i] >= MX1 and X1[i] <= max(X1):
           Xn.insert(j, X1[i])
           Yn.insert(j, Y1[i])
           j=j+1

    for i in range(0, len(X1)):
        if X1[i] > mX1 and X1[i] < MX1 and Y1[i] < 0 and abs(Y1[i]-med_do) < 1.0*Sigma_do:
           Yn.insert(j, Y1[i])
           Xn.insert(j, X1[i])
           j=j+1

    Xn=np.array(Xn)
    Yn=np.array(Yn)

    # Rotation and translation to put in non-canonical system
    X2 = x0 + Xn * np.cos(phi) - Yn * np.sin(phi)
    Y2 = y0 + Xn * np.sin(phi) + Yn * np.cos(phi) 

    return X2, Y2

#==================================================================================

def rms_superellipse_pts(params, X, Y):

    """
    Compute rms of the points X, Y using the superellipse with n=6 
    described by the params = x0, y0, ap, bp, phi.
    A. Carbognani, INAF-OAS, July 27, 2022
    """

    x0, y0, ap, bp, phi = params

    # Canonical superellipse reduction
    X1=((X-x0)*np.cos(phi)+(Y-y0)*np.sin(phi))/ap;
    Y1=((Y-y0)*np.cos(phi)-(X-x0)*np.sin(phi))/bp;
    
    # Algebric distance between observed points and superellipse equation with n=6
    Q=(X1*X1*X1*X1*X1*X1 + Y1*Y1*Y1*Y1*Y1*Y1 - 1); 

    # Compute square of distances
    Q_sqr=Q*Q;

    # Compute mean rms for points 
    rms=np.sqrt(sum(Q_sqr))/len(Q_sqr)

    return rms

#==================================================================================

def track_center_superellipse(x, y): 
    
    """
    Compute best track's center using the best fit superellipse with n=6 
    described by the params = x0, y0, ap, bp, phi.
    A. Carbognani, INAF-OAS, July 27, 2022
    """

    # Compute best fit ellipse with raw data
    coeffs = fit_ellipse(x, y)
    x0, y0, ap, bp, e, phi = cart_to_pol(coeffs)
        
    # Smooth raw data  
    x2, y2 = smooth_superellipse_pts((ap, bp, phi), x, y) 
    
    # Compute best fit ellipse with smooth data
    coeffs1 = fit_ellipse(np.transpose(x2), np.transpose(y2))
    x0, y0, ap, bp, e, phi = cart_to_pol(coeffs1)

    # Smooth already smooth data using second guest superellipse with n=4
    #x2, y2 = smooth_superellipse_pts((ap, bp, phi), x1, y1) 
    
    # Compute best fit ellipse with smooth-smooth data
    #coeffs2 = fit_ellipse(np.transpose(x2), np.transpose(y2))
    #x0, y0, ap, bp, e, phi = cart_to_pol(coeffs2)

    # Plot guest superellipse  
    x1, y1 = get_superellipse_pts((x0, y0, ap, bp, phi)) 
    x22, y22 = get_ellipse_pts((x0, y0, ap, bp, e, phi))
    plt.plot(x1, y1, label = 'Guest superellipse')
    plt.plot(x22, y22, label = 'Best fit ellipse')

    #################################
    # Compute best fit superellipse #
    #################################  

    # Compute rms for guest-superellipse
    rmsold=rms_superellipse_pts((x0, y0, ap, bp, phi), x2, y2)
    print('First rms', rmsold)

    # Set grid resolution
    ResGrid=3;
    
    # Starting domains of variation of the superellipse's parameters (phi is fixed)
    var=0.2
    dXdom=var*x0
    dYdom=var*y0
    dapdom=var*ap
    dbpdom=var*bp
    
    # Minimum rms array index
    minrms=[0, 0, 0, 0]
    
    while dXdom > 0.5:
        for i in range(-ResGrid, ResGrid+1):
            xx0=x0+i*dXdom/ResGrid
            for j in range(-ResGrid, ResGrid+1):  
                yy0=y0+j*dYdom/ResGrid
                for k in range(-ResGrid, ResGrid+1): 
                    aap=ap+k*dapdom/ResGrid
                    for l in range(-ResGrid, ResGrid+1):
                        bbp=bp+l*dbpdom/ResGrid
                        if aap > bbp:
                            rms=rms_superellipse_pts((xx0, yy0, aap, bbp, phi), x2, y2)
                            if rms<rmsold:
                                minrms.insert(0, i) # Save best index for x0
                                minrms.insert(1, j) # Save best index for y0
                                minrms.insert(2, k) # Save best index for ap
                                minrms.insert(3, l) # Save best index for bp
                                rmsold=rms

        print('second', rmsold)
        # New best fit superellipse parameters with n=4
        x0=x0+minrms[0]*dXdom/ResGrid
        y0=y0+minrms[1]*dYdom/ResGrid
        ap=ap+minrms[2]*dapdom/ResGrid
        bp=bp+minrms[3]*dbpdom/ResGrid
           
        # Dimezzamento domini di ricerca
        dXdom=dXdom/2
        dYdom=dYdom/2
        dapdom=dapdom/2
        dbpdom=dbpdom/2

    # Final plots
    # Plot original points
    plt.plot(x, y, 'x', label = 'Original points')
    x, y = get_superellipse_pts((x0, y0, ap, bp, phi)) # Best-fit superellipse
    plt.plot(x, y, label = 'Best fit superellipse')
    plt.plot(x2, y2, 'o', label = 'Smooth points')   # smooth points
    # plt.gca().set_aspect('equal', adjustable='box') # Equal axis scale
    plt.legend(loc="upper left")
    plt.show()
    
    # Best fit track's center
    return x0, y0

#==================================================================================

def track_center(x, y, phi):
 
    """
    Compute best track's center using border's median algorithm
    x, y are contours array, phi is the inclination angle of the track (radiant)
    A. Carbognani, INAF-OAS, Aug 03, 2022
    """

    # Put the best fit track's center 
    # to make them less sensitive to edge irregularities.
    x0=statistics.mean(x)
    y0=statistics.mean(y)

    # Compute best track's center with border's median algorithm
    x1, y1, X0, Y0 = smooth_borders_median((x0, y0, phi), x, y) 
    
    #print('Median center', X0, Y0)

    # Plot original points
    #plt.plot(x, y, 'x', label = 'Original points')
    #plt.plot(x1, y1, '.', label = 'Best fit track')   # Median borders
    #plt.legend(loc="upper left")
    # plt.gca().set_aspect('equal', adjustable='box') # Equal axis scale
    #plt.show()

    # Best fit track's center
    return X0, Y0
#==================================================================================

def track_center2(x, y):
 
    """
    Compute best track's center first deleting outliers and after using border's median algorithm
    x, y are raw contours array
    A. Carbognani, INAF-OAS, Aug 04, 2022
    """

    # Compute best fit ellipse with track's raw data
    coeffs = fit_ellipse(x, y)

    # Geometrical parameter best fit ellipse
    x0, y0, ap, bp, e, phi = cart_to_pol(coeffs)

    # Compute best-fit ellipse's perimeter 
    x2, y2 = get_ellipse_pts((x0, y0, ap, bp, e, phi))        

    # Delete outliers
    x1, y1 = smooth_borders_median2((x0, y0, phi), x, y) 

    # Compute best fit ellipse with no-outliers data
    coeffs = fit_ellipse(x1, y1)
    
    # Geometrical parameter second best fit ellipse
    x0, y0, ap, bp, e, phi = cart_to_pol(coeffs)

    # Compute second best-fit ellipse's perimeter 
    x2, y2 = get_ellipse_pts((x0, y0, ap, bp, e, phi))  

    # Compute best fit track's center with border median algorithm
    x3, y3, X0, Y0 = smooth_borders_median((x0, y0, phi), x1, y1) 

    # Plot original points
    #plt.plot(x, y, 'x', label = 'Original points')
    # Plot best-fit ellipse
    #plt.plot(x2, y2, label = 'Best fit ellipse with no-outliers')
    #plt.plot(x1, y1, '.', label = 'No outliers border')  
    #plt.plot(x3, y3, '.', label = 'Best fit track')   # Median borders
    #plt.legend(loc="upper left")
    # plt.gca().set_aspect('equal', adjustable='box') # Equal axis scale
    #plt.show()

    # Best fit track's center
    return X0, Y0

#==================================================================================

def track_center_median(x, y):
 
    """
    Compute best track's center using median function
    A. Carbognani, INAF-OAS, July 29, 2022
    """

    x0=statistics.median(x)
    y0=statistics.median(y)

    # Best fit track's center
    return x0, y0

#==================================================================================

def track_center_mean(x, y):
 
    """
    Compute best track's center using mean function
    A. Carbognani, INAF-OAS, July 29, 2022
    """

    x0=statistics.mean(x)
    y0=statistics.mean(y)

    # Best fit track's center
    return x0, y0

#==================================================================================



















