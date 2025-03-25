from tkinter import *
from tkinter import ttk
import sv_ttk
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk) 
from matplotlib import *
from numpy import *
import math
import numpy as np 
import simpy as sp
import scipy.spatial as spatial
import scipy.integrate as integrate
import pathlib as path

# list of squares
def vortexpanel(x,y):
    class Pannello:
    
        def __init__(self, xa, ya, xb, yb):
        
            self.xa, self.ya = xa, ya  # start
            self.xb, self.yb = xb, yb  # end
        
            self.xc, self.yc = (xa + xb) / 2, (ya + yb) / 2  #  center
            self.length = np.sqrt((xb - xa)**2 + (yb - ya)**2)  # panel length
            
            # orientation
            if xb - xa <= 0.0:
                self.beta = np.arccos((yb - ya) / self.length)
            elif xb - xa > 0.0:
                self.beta = np.pi + np.arccos(-(yb - ya) / self.length)
            
            # panel location
            if self.beta <= np.pi:
                self.loc = 'upper'  # upper surface
            else:
                self.loc = 'lower'  # lower surface
            
            self.sigma = 0.0  # source strength
            self.vt = 0.0  # tangential velocity
            self.cp = 0.0  # pressure coefficient
    def define_panels(x, y, N=40):
        R = (x.max() - x.min()) / 2.0  # circle radius
        x_center = (x.max() + x.min()) / 2.0  # x-coordinate of circle center
        
        teta = np.linspace(0.0, 2.0 * np.pi, N + 1)  # array of angles
        x_circle = x_center + R * np.cos(teta)  # x-coordinates of circle
        
        x_ends = np.copy(x_circle)  # x-coordinate of panels end-points
        y_ends = np.empty_like(x_ends)  # y-coordinate of panels end-points
        
        # extend coordinates to consider closed surface
        x, y = np.append(x, x[0]), np.append(y, y[0])
        
        # compute y-coordinate of end-points by projection
        I = 0
        for i in range(N):
            while I < len(x) - 1:
                if (x[I] <= x_ends[i] <= x[I + 1]) or (x[I + 1] <= x_ends[i] <= x[I]):
                    break
                else:
                    I += 1
            a = (y[I + 1] - y[I]) / (x[I + 1] - x[I])
            b = y[I + 1] - a * x[I + 1]
            y_ends[i] = a * x_ends[i] + b
        y_ends[N] = y_ends[0]
        
        # create panels
        panels = np.empty(N, dtype=object)
        for i in range(N):
            panels[i] = Pannello(x_ends[i], y_ends[i], x_ends[i + 1], y_ends[i + 1])
        
        return panels
    panels = define_panels(x, y, N=40)
    class Freestream:
        def __init__(self, u_inf=1.0, alpha=0.0):
            
            self.u_inf = u_inf
            self.alpha = np.radians(alpha)
    #poi cambia per mettere i tuoi angoli di attacco e di forza
    freestream = Freestream(u_inf=1.0, alpha=4.0)
    def integral(x, y, panel, dxdk, dydk):

        def integrand(s):
            return (((x - (panel.xa - np.sin(panel.beta) * s)) * dxdk +
                    (y - (panel.ya + np.cos(panel.beta) * s)) * dydk) /
                    ((x - (panel.xa - np.sin(panel.beta) * s))**2 +
                    (y - (panel.ya + np.cos(panel.beta) * s))**2) )
        
        return integrate.quad(integrand, 0.0, panel.length)[0]
    def source_contribution_normal(panels):
        # crea la prima matrice di risoluzione 
        A = np.empty((panels.size, panels.size), dtype=float)
        # source contribution on a panel from itself
        np.fill_diagonal(A, 0.5)
        # source contribution on a panel from others
        for i, panel_i in enumerate(panels):
            for j, panel_j in enumerate(panels):
                if i != j:
                    A[i, j] = 0.5 / np.pi * integral(panel_i.xc, panel_i.yc, 
                                                        panel_j,
                                                        np.cos(panel_i.beta),
                                                        np.sin(panel_i.beta))
        return A
    def vortex_contribution_normal(panels):
   
        A = np.empty((panels.size, panels.size), dtype=float)
        # vortex contribution on a panel from itself
        np.fill_diagonal(A, 0.0)
        # vortex contribution on a panel from others
        for i, panel_i in enumerate(panels):
            for j, panel_j in enumerate(panels):
                if i != j:
                    A[i, j] = -0.5 / np.pi * integral(panel_i.xc, panel_i.yc, 
                                                        panel_j,
                                                        np.sin(panel_i.beta),
                                                        -np.cos(panel_i.beta))
        return A
    A_source = source_contribution_normal(panels)
    B_vortex = vortex_contribution_normal(panels)

    def kutta_condition(A_source, B_vortex):

    #Builds the Kutta condition array.
    
 
        b = np.empty(A_source.shape[0] + 1, dtype=float)
        # matrix of source contribution on tangential velocity
        # is the same than
        # matrix of vortex contribution on normal velocity
        b[:-1] = B_vortex[0, :] + B_vortex[-1, :]
        # matrix of vortex contribution on tangential velocity
        # is the opposite of
        # matrix of source contribution on normal velocity
        b[-1] = - np.sum(A_source[0, :] + A_source[-1, :])
        return b
    def build_singularity_matrix(A_source, B_vortex):

        A = np.empty((A_source.shape[0] + 1, A_source.shape[1] + 1), dtype=float)
        # source contribution matrix
        A[:-1, :-1] = A_source
        # vortex contribution array
        A[:-1, -1] = np.sum(B_vortex, axis=1)
        # Kutta condition array
        A[-1, :] = kutta_condition(A_source, B_vortex)
        return A
    def build_freestream_rhs(panels, freestream):
   
        b = np.empty(panels.size + 1, dtype=float)
        # freestream contribution on each panel
        for i, panel in enumerate(panels):
            b[i] = -freestream.u_inf * np.cos(freestream.alpha - panel.beta)
        # freestream contribution on the Kutta condition
        b[-1] = -freestream.u_inf * (np.sin(freestream.alpha - panels[0].beta) +
                                    np.sin(freestream.alpha - panels[-1].beta) )
        return b
    A = build_singularity_matrix(A_source, B_vortex)
    b = build_freestream_rhs(panels, freestream)
    strengths = np.linalg.solve(A, b)
    # store source strength on each panel
    for i , panel in enumerate(panels):
        panel.sigma = strengths[i]
        
    # store circulation density
    gamma = strengths[-1]
    def compute_tangential_velocity(panels, freestream, gamma, A_source, B_vortex):
        A = np.empty((panels.size, panels.size + 1), dtype=float)
        # matrix of source contribution on tangential velocity
        # is the same than
        # matrix of vortex contribution on normal velocity
        A[:, :-1] = B_vortex
        # matrix of vortex contribution on tangential velocity
        # is the opposite of
        # matrix of source contribution on normal velocity
        A[:, -1] = -np.sum(A_source, axis=1)
        # freestream contribution
        b = freestream.u_inf * np.sin([freestream.alpha - panel.beta 
                                        for panel in panels])
        
        strengths = np.append([panel.sigma for panel in panels], gamma)
        
        tangential_velocities = np.dot(A, strengths) + b
        
        for i, panel in enumerate(panels):
            panel.vt = tangential_velocities[i]
    compute_tangential_velocity(panels, freestream, gamma, A_source, B_vortex)
    def compute_pressure_coefficient(panels, freestream):
        for panel in panels:
            panel.cp = 1.0 - (panel.vt / freestream.u_inf)**2
    compute_pressure_coefficient(panels, freestream)
    c = abs(max(panel.xa for panel in panels) - min(panel.xa for panel in panels))
    cl = (gamma * sum(panel.length for panel in panels) /  (0.5 * freestream.u_inf * c))

def STREAMLINE_VPM(XP,YP,XB,YB,phi,S):
      # Number of panels
    numPan = len(XB)-1                                                          # Number of panels (control points)
    
    # Initialize arrays
    Nx = np.zeros(numPan)                                                       # Initialize Nx integral array
    Ny = np.zeros(numPan)                                                       # Initialize Ny integral array
    
    # Compute Nx and Ny
    for j in range(numPan):                                                     # Loop over all panels
        # Compute intermediate values
        A = -(XP-XB[j])*np.cos(phi[j]) - (YP-YB[j])*np.sin(phi[j])              # A term
        B  = (XP-XB[j])**2 + (YP-YB[j])**2                                      # B term
        Cx = np.sin(phi[j])                                                     # Cx term (X-direction)
        Dx = -(YP-YB[j])                                                        # Dx term (X-direction)
        Cy = -np.cos(phi[j])                                                    # Cy term (Y-direction)
        Dy = XP-XB[j]                                                           # Dy term (Y-direction)
        E  = math.sqrt(B-A**2)                                                  # E term
        if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):           # If E term is 0 or complex or a NAN or an INF
            Nx[j] = 0                                                           # Set Nx value equal to zero
            Ny[j] = 0                                                           # Set Ny value equal to zero
        else:
            # Compute Nx, Ref [1]
            term1 = 0.5*Cx*np.log((S[j]**2 + 2*A*S[j]+B)/B);                    # First term in Nx equation
            term2 = ((Dx-A*Cx)/E)*(math.atan2((S[j]+A),E) - math.atan2(A,E));   # Second term in Nx equation
            Nx[j] = term1 + term2;                                              # Compute Nx integral
            
            # Compute Ny, Ref [1]
            term1 = 0.5*Cy*np.log((S[j]**2 + 2*A*S[j]+B)/B);                    # First term in Ny equation
            term2 = ((Dy-A*Cy)/E)*(math.atan2((S[j]+A),E) - math.atan2(A,E));   # Second term in Ny equation
            Ny[j] = term1 + term2;                                              # Compute Ny integral
            
        # Zero out any problem values
        if (np.iscomplex(Nx[j]) or np.isnan(Nx[j]) or np.isinf(Nx[j])):         # If Nx term is complex or a NAN or an INF
            Nx[j] = 0                                                           # Set Nx value equal to zero
        if (np.iscomplex(Ny[j]) or np.isnan(Ny[j]) or np.isinf(Ny[j])):         # If Ny term is complex or a NAN or an INF
            Ny[j] = 0                                                           # Set Ny value equal to zero
    
    return Nx, Ny                                                               # Return both Nx and Ny matrices

def streamlines():
        # Streamline parameters like source panels, ma poi aggiungi la roba da streamline vp
    nGridX   = 150                                                              # X-grid for streamlines and contours
    nGridY   = 150                                                              # Y-grid for streamlines and contours
    xVals    = [-0.5, 1.5]                                                      # X-grid extents [min, max]
    yVals    = [-0.3, 0.3]     
                                                     
    Xgrid  = np.linspace(xVals[0],xVals[1],nGridX)                              # X-values in evenly spaced grid
    Ygrid  = np.linspace(yVals[0],yVals[1],nGridY)                              # Y-values in evenly spaced grid
    XX, YY = np.meshgrid(Xgrid,Ygrid)   

    Ysl    = np.linspace(yVals[0],yVals[1],int(0.25*nGridY))             # Create array of Y streamline starting points
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
    CpXY = 1 - (Vxy/Vinf)**2    

def plot1(): 
    
    # the figure that will contain the plot 
    fig = Figure(figsize = (5, 5), 
                 dpi = 100) 
  
    # adding the subplot 
    plot1 = fig.add_subplot(111) 
    plot1.set_xlim(0,100)
    #impostiamo i dati presi dalla GUI, e creaiamo un linspace
    corda= 1
    x=np.array(linspace(0, float(corda), num=1000))
    naca=nacaa_var.get()
    #dividiamo naca cifra per cifra
    arr_naca=list(naca)
    naca_subtype = len(arr_naca)
    #normalizziamo per la corda
    xn=x/float(corda)
    #creaiamo 3 casi: naca 4 cifre simmetrichi, naca 4 cifre non simmetrici, naca 5 cifre
    if naca_subtype == 4:
        #impostiamo i valori con cui trovare il profilo
        m=(int(arr_naca[0]))/100
        p=int(arr_naca[1])/10
        t=(10*int(arr_naca[2])+int(arr_naca[3]))/100
        #troviamo l'equazione della mediana
        if m==0 and p==0:
            # troviamo equazione della parte superiore e della parte inferiore del profilo
            yt=(t/0.2)*((0.2969*xn**(0.5))-(0.126*xn)-(0.3516*xn**2)+(0.2843*xn**3)-(0.1015*xn**4))
            xu=np.array(xn)
            xl=np.array(xn)
            yu=np.array(yt)
            yl=np.array(-yt)
            ymediana=0*x
            xf=np.concatenate((np.flip(xu), xl))
            yf=np.concatenate((np.flip(yu), yl))
        else:
            # troviamo equazione della parte superiore e della parte inferiore del profilo
            yt=(t/0.2)*((0.2969*xn**(0.5))-(0.126*xn)-(0.3516*xn**2)+(0.2843*xn**3)-(0.1015*xn**4))
            for i in xn:
                #troviamo l'equazione della mediana
                if i<p:
                    ymediana=(m/(p**2))*(2*p*xn-xn**2)
                elif i>p:
                    ymediana=(m/(1-p)**2)*((1-2*p)+2*p*xn-xn**2)

            # troviamo equazione della parte superiore e della parte inferiore del profilo
            teth=np.gradient(ymediana,xn)
            xu=np.array(xn-yt*sin(teth))
            xl=np.array(xn+yt*sin(teth))
            yu=np.array(ymediana+yt*cos(teth))
            yl=np.array(ymediana-yt*cos(teth))
            xf=np.concatenate((np.flip(xu), xl))
            yf=np.concatenate((np.flip(yu), yl))
    elif naca_subtype == 5:
        cl=1.5*int(arr_naca[0])
        if int(arr_naca[1])<=0:
            return None
        if int(arr_naca[1])>6:
            return None
        
        p=0.5*(10*int(arr_naca[1])+int(arr_naca[2]))
        t=10*int(arr_naca[3])+int(arr_naca[4])
        m_arr=np.array([0.058,0.126,0.2025,0.29,0.391])
        k_one_arr=np.array([361.4,51.64,15.957,6.643,3.23])
        arr_naca_int=int(arr_naca[1])
        m=float(m_arr[arr_naca_int-1])
        k_one=float(k_one_arr[arr_naca_int-1])
        yt=(t/0.2)*((0.2969*xn**(0.5))-(0.126*xn)-(0.3516*xn**2)+(0.2843*xn**3)-(0.1015*xn**4))

        for i in xn:
            if i<p:
                ymediana=(k_one/6)*(xn**3-3*m*xn**2+xn*(3-m)*m**2)
            elif i>p:
                ymediana=(k_one*(m**3)/6)*(1-xn) 
        teth=np.gradient(ymediana,xn)
        xu=np.array(xn-yt*sin(teth))
        xl=np.array(xn+yt*sin(teth))
        yu=np.array(ymediana+yt*cos(teth))
        yl=np.array(ymediana-yt*cos(teth))
        xf=np.concatenate((np.flip(xu), xl))
        yf=np.concatenate((np.flip(yu), yl))
    else:
        return None  # Restituisce None se il numero non ha 4 o 5 arr_naca
    # definiamo spazio del plot
    vortexpanels(xf,yf)
    plot1.set_xlim(-1.5,1.5)
    plot1.set_ylim(-max(xu)*1.5,max(xu)*1.5)
    plot1.plot(xu,yu,'k-',xl,yl,'k-',x,ymediana,'r-',x_ends,y_ends) 

    # creating the Tkinter canvas 
    # containing the Matplotlib figure 
    canvas = FigureCanvasTkAgg(fig, 
                               master = window)   
    canvas.get_tk_widget().grid(row=5,column=1,columnspan=3,rowspan=20,sticky=N)
    # here: plot suff to your fig
    canvas.draw()
  
def calculate(*args):
    try:
        value = float(nacaa_var.get())
        chord_value.set(float(value*0.1))
    except ValueError:
        pass

# the main Tkinter window 
window = Tk() 
w = Scrollbar(window) 
window.columnconfigure(3)
window.rowconfigure(3)
# setting the title 
window.title('Plotting in Tkinter') 

# dimensions of the main window 

# buttons that displays the 3 plots
plot_button = Button(master = window, 
					command = plot1, 
					height = 10, 
					width = 20, 
					text = "Jukowsky") 
plot_button2 = Button(master = window, 
					command = plot1, 
					height = 10, 
					width = 20, 
					text = "Profili Sottili") 
plot_button3 = Button(master = window, 
					command = plot1, 
					height = 10, 
					width = 20, 
					text = "Metodo a Pannelli") 

# place the buttons 
# in main window 
plot_button.grid(row=1,column=1,sticky=N)
plot_button2.grid(row=1,column=2,sticky=N)
plot_button3.grid(row=1,column=3,sticky=N) 

#piazzare input di testo in cui mettere i dati
nacaa_var=StringVar()
nacaa_entry = ttk.Entry(window, width=7, textvariable=nacaa_var)
nacaa_entry.grid(column=2, row=2, sticky=(W, E))

chord_value=StringVar()
curv_entry = ttk.Entry(window, width=7, textvariable=chord_value)
curv_entry.grid(column=3, row=2, sticky=(W, E))

nacaa_entry.focus()
window.bind("<Return>", calculate())

#stile un po' più estetico
sv_ttk.use_dark_theme()
# run the gui 
window.mainloop()