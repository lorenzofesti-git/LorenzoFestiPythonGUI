import numpy as np
import math as math
import matplotlib.pyplot as plt
from matplotlib import path
from import2 import COMPUTE_KL_VPM
from import3 import STREAMLINE_VPM
from import4 import COMPUTE_CIRCULATION

def pannelli(Vinf,AoA,XF,YF):
    # Vinf = 1                                                                        # Freestream velocity [] (just leave this at 1)
    # AoA  = 0                                                                        # Angle of attack [deg]
    # Convert angle of attack to radians
    AoAR = AoA*(np.pi/180)      
    for i in range(0,200):
    # inti=int(i)*10
        XB[i]=XF[i*10]
        YB[i]=YF[i*10]


    numPts = len(XB)                                                                # Number of boundary points
    numPan = numPts - 1                                                             # Number of panels (control points)

    edge = np.zeros(numPan)                                                         # Initialize edge value array
    for i in range(numPan):                                                         # Loop over all panels
        edge[i] = (XB[i+1]-XB[i])*(YB[i+1]+YB[i])                                   # Compute edge values

    sumEdge = np.sum(edge)                                                          # Sum all edge values

    # If panels are CCW, flip them (don't if CW)
    if (sumEdge < 0):                                                               # If panels are CCW
        XB = np.flipud(XB)                                                          # Flip the X-data array
        YB = np.flipud(YB)                                                          # Flip the Y-data array

    # Initialize variables
    XC  = np.zeros(numPan)                                                          # Initialize control point X-coordinate array
    YC  = np.zeros(numPan)                                                          # Initialize control point Y-coordinate array
    S   = np.zeros(numPan)                                                          # Initialize panel length array
    phi = np.zeros(numPan)                                                          # Initialize panel orientation angle array [deg]

    # Find geometric quantities of the airfoil
    for i in range(numPan):                                                         # Loop over all panels
        XC[i]   = 0.5*(XB[i]+XB[i+1])                                               # X-value of control point
        YC[i]   = 0.5*(YB[i]+YB[i+1])                                               # Y-value of control point
        dx      = XB[i+1]-XB[i]                                                     # Change in X between boundary points
        dy      = YB[i+1]-YB[i]                                                     # Change in Y between boundary points
        S[i]    = (dx**2 + dy**2)**0.5                                              # Length of the panel
        phi[i]  = math.atan2(dy,dx)                                                 # Angle of panel (positive X-axis to inside face)
        if (phi[i] < 0):                                                            # Make all panel angles positive [rad]
            phi[i] = phi[i] + 2*np.pi

    # Compute angle of panel normal w.r.t. horizontal and include AoA
    delta                = phi + (np.pi/2)                                          # Angle from positive X-axis to outward normal vector [rad]
    beta                 = delta - AoAR                                             # Angle between freestream vector and outward normal vector [rad]
    beta[beta > 2*np.pi] = beta[beta > 2*np.pi] - 2*np.pi                           # Make all panel angles between 0 and 2pi [rad]

    # Geometric integral (normal [K] and tangential [L])
    # - Refs [2] and [3]
    K, L = COMPUTE_KL_VPM(XC,YC,XB,YB,phi,S)                                        # Compute geometric integrals

    # Populate A matrix
    A = np.zeros([numPan,numPan])                                                   # Initialize the A matrix
    for i in range(numPan):                                                         # Loop over all i panels
        for j in range(numPan):                                                     # Loop over all j panels
            if (i == j):                                                            # If the panels are the same
                A[i,j] = 0                                                          # Set A equal to zero
            else:                                                                   # If panels are not the same
                A[i,j] = -K[i,j]                                                    # Set A equal to negative geometric integral
                
    # Populate b array
    b = np.zeros(numPan)                                                            # Initialize the b array
    for i in range(numPan):                                                         # Loop over all panels
        b[i] = -Vinf*2*np.pi*np.cos(beta[i])                                        # Compute RHS array

    # Satisfy the Kutta condition
    pct    = 100                                                                    # Panel replacement percentage
    panRep = int((pct/100)*numPan)-1                                                # Replace this panel with Kutta condition equation
    if (panRep >= numPan):                                                          # If we specify the last panel
        panRep = numPan-1                                                           # Set appropriate replacement panel index
    A[panRep,:]        = 0                                                          # Set all colums of the replaced panel equation to zero
    A[panRep,0]        = 1                                                          # Set first column of replaced panel equal to 1
    A[panRep,numPan-1] = 1                                                          # Set last column of replaced panel equal to 1
    b[panRep]          = 0                                                          # Set replaced panel value in b array equal to zero

    # Compute gamma values
    gamma = np.linalg.solve(A,b)                                                    # Compute all vortex strength values

    # Compute velocities
    Vt = np.zeros(numPan)                                                           # Initialize tangential velocity array
    Cp = np.zeros(numPan)                                                           # Initialize pressure coefficient array
    for i in range(numPan):                                                         # Loop over all i panels
        addVal = 0                                                                  # Reset summation value to zero
        for j in range(numPan):                                                     # Loop over all j panels
            addVal = addVal - (gamma[j]/(2*np.pi))*L[i,j]                           # Sum all tangential vortex panel terms
        
        Vt[i] = Vinf*np.sin(beta[i]) + addVal + gamma[i]/2                          # Compute tangential velocity by adding uniform flow and i=j terms
        Cp[i] = 1 - (Vt[i]/Vinf)**2                                                 # Compute pressure coefficient

    # Compute normal and axial force coefficients
    CN = -Cp*S*np.sin(beta)                                                         # Normal force coefficient []
    CA = -Cp*S*np.cos(beta)                                                         # Axial force coefficient []

    # Compute lift, drag, and moment coefficients
    CL = sum(CN*np.cos(AoAR)) - sum(CA*np.sin(AoAR))                                # Decompose axial and normal to lift coefficient []
    CM = sum(Cp*(XC-0.25)*S*np.cos(phi))                                            # Moment coefficient []
                                    # If we are plotting streamlines or pressure coefficient contours
    # Grid parameters
    nGridX   = 150                                                              # X-grid for streamlines and contours
    nGridY   = 150                                                              # Y-grid for streamlines and contours
    xVals    = [-0.5, 1.5]                                                      # X-grid extents [min, max]
    yVals    = [-0.3, 0.3]                                                      # Y-grid extents [min, max]

    # Streamline parameters
    slPct  = 25                                                                 # Percentage of streamlines of the grid
    Ysl    = np.linspace(yVals[0],yVals[1],int((slPct/100)*nGridY))             # Create array of Y streamline starting points
    Xsl    = xVals[0]*np.ones(len(Ysl))                                         # Create array of X streamline starting points
    XYsl   = np.vstack((Xsl.T,Ysl.T)).T                                         # Concatenate X and Y streamline starting points

    # Generate the grid points
    Xgrid  = np.linspace(xVals[0],xVals[1],nGridX)                              # X-values in evenly spaced grid
    Ygrid  = np.linspace(yVals[0],yVals[1],nGridY)                              # Y-values in evenly spaced grid
    XX, YY = np.meshgrid(Xgrid,Ygrid)                                           # Create meshgrid from X and Y grid arrays

    # Initialize velocities
    Vx     = np.zeros([nGridX,nGridY])                                          # Initialize X velocity matrix
    Vy     = np.zeros([nGridX,nGridY])                                          # Initialize Y velocity matrix

    # Path to figure out if grid point is inside polygon or not
    AF     = np.vstack((XB.T,YB.T)).T                                           # Concatenate XB and YB geometry points
    afPath = path.Path(AF)                                                      # Create a path for the geometry

    # Solve for grid point X and Y velocities
    for m in range(nGridX):                                                     # Loop over X-grid points
        for n in range(nGridY):                                                 # Loop over Y-grid points
            XP     = XX[m,n]                                                    # Current iteration's X grid point
            YP     = YY[m,n]                                                    # Current iteration's Y grid point
            Nx, Ny = STREAMLINE_VPM(XP,YP,XB,YB,phi,S)                          # Compute Nx and Ny geometric integrals
            
            # Check if grid points are in object
            # - If they are, assign a velocity of zero
            if afPath.contains_points([(XP,YP)]):                               # If (XP,YP) is in the body
                Vx[m,n] = 0                                                     # Set X-velocity equal to zero
                Vy[m,n] = 0                                                     # Set Y-velocity equal to zero
            else:
                Vx[m,n] = Vinf*np.cos(AoAR) + sum(-gamma*Nx/(2*np.pi))          # Compute X-velocity
                Vy[m,n] = Vinf*np.sin(AoAR) + sum(-gamma*Ny/(2*np.pi))          # Compute Y-velocity

    # Compute grid point velocity magnitude and pressure coefficient
    Vxy  = np.sqrt(Vx**2 + Vy**2)                                               # Compute magnitude of velocity vector []
    CpXY = 1 - (Vxy/Vinf)**2                                                    # Pressure coefficient []


     # If we are plotting streamlines or Cp contours
    # Compute circulation
    aa   = 0.75                                                                 # Ellipse horizontal half-length
    bb   = 0.25                                                                 # Ellipse vertical half-length
    x0   = 0.5                                                                  # Ellipse center X-coordinate
    y0   = 0                                                                    # Ellipse center Y-coordinate
    numT = 5000                                                                 # Number of points on ellipse
    Circulation, xC, yC, VxC, VyC = COMPUTE_CIRCULATION(aa,bb,x0,y0,            # Compute circulation around ellipse
                                                        numT,Vx,Vy,Xgrid,Ygrid)
