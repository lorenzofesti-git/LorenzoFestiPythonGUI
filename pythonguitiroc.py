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
from panelmethod import *
from thinairf_solve import *
from Joukoski_transf_curves import *
from nacanumber import *

def met_panels(): 
    #impostiamo i dati presi dalla GUI, e creaiamo un linspace
        # the figure that will contain the plot 
    fig = Figure(figsize = (5,4), dpi = 100)  
    # adding the subplot 
    met_panels = fig.add_subplot(111) 
    corda= 1
    x=np.array(linspace(0, float(corda), num=100))
    naca=nacaa_var.get()

    xu,xl,yu,yl,xf,yf,ymediana,deltY= nacanum(naca,x)
    
    xb,yb,VxC, VyC=pannelli(float(U_inf_value),float(AoA_value),xf,yf)


    # definiamo spazio del plot

        # creating the Tkinter canvas 
    # containing the Matplotlib figure 
    canvas = FigureCanvasTkAgg(fig, 
                                master = window)   
    canvas.get_tk_widget().grid(row=5,column=1,columnspan=3,rowspan=20,sticky=N)
        # here: plot suff to your fig

    canvas.draw()
    
    met_panels.set_xlim(-0.5,1.2)
    met_panels.set_ylim(-max(xu)*1.2,max(xu)*1.2)
    met_panels.plot(xu,yu,'k-',xl,yl,'k-',x,ymediana,'r-',xb,yb,'2') 
    #,xb,yb,'2'

def met_prof_sott(): 
    #impostiamo i dati presi dalla GUI, e creaiamo un linspace
        # the figure that will contain the plot 
    fig = Figure(figsize = (5,4), dpi = 100)  
    # adding the subplot 
    met_panels = fig.add_subplot(111) 
    corda= 1
    x=np.array(linspace(0, float(corda), num=100))
    naca=nacaa_var.get()
    # xb=zeros(shape=(200,1))
    # yb=zeros(shape=(200,1))
    xu,xl,yu,yl,xf,yf,ymediana,deltY= nacanum(naca,x)

    A0, A1, A2 = solve_thin_airf_theory(angle_of_attack=5*np.pi/180, n_coefficients=3)


    # definiamo spazio del plot

        # creating the Tkinter canvas 
    # containing the Matplotlib figure 
    canvas = FigureCanvasTkAgg(fig, 
                                master = window)   
    canvas.get_tk_widget().grid(row=5,column=1,columnspan=3,rowspan=20,sticky=N)
        # here: plot suff to your fig

    canvas.draw()
    
    met_panels.set_xlim(-0.5,1.2)
    met_panels.set_ylim(-max(xu)*1.2,max(xu)*1.2)
    met_panels.plot(xu,yu,'k-',xl,yl,'k-',x,ymediana,'r-',xb,yb,'2') 
    #,xb,yb,'2'
#

def met_trasf_conf(): 
    #impostiamo i dati presi dalla GUI, e creaiamo un linspace
        # the figure that will contain the plot 
    fig = Figure(figsize = (5,4), dpi = 100)  
    # adding the subplot 
    met_panels = fig.add_subplot(111) 
    corda= 1
    x=np.array(linspace(0, float(corda), num=100))
    naca=nacaa_var.get()
    # xb=zeros(shape=(200,1))
    # yb=zeros(shape=(200,1))
    xu,xl,yu,yl,xf,yf,ymediana,deltY= nacanum(naca,x)

    A0, A1, A2 = solve_thin_airf_theory(angle_of_attack=5*np.pi/180, n_coefficients=3)


    # definiamo spazio del plot

        # creating the Tkinter canvas 
    # containing the Matplotlib figure 
    canvas = FigureCanvasTkAgg(fig, 
                                master = window)   
    canvas.get_tk_widget().grid(row=5,column=1,columnspan=3,rowspan=20,sticky=N)
        # here: plot suff to your fig

    canvas.draw()
    
    met_panels.set_xlim(-0.5,1.2)
    met_panels.set_ylim(-max(xu)*1.2,max(xu)*1.2)
    met_panels.plot(xu,yu,'k-',xl,yl,'k-',x,ymediana,'r-',xb,yb,'2') 
    #,xb,yb,'2'


# the main Tkinter window 
window = Tk() 
window.columnconfigure(3)
window.rowconfigure(3)

# setting the title 
window.title('Plotting') 

# dimensions of the main window 

# buttons that displays the 3 plots
plot_button = Button(master = window, 
					command = met_panels, 
					height = 10, 
					width = 20, 
					text = "Metodo a Pannelli") 
plot_button2 = Button(master = window, 
					command = met_prof_sott, 
					height = 10, 
					width = 20, 
					text = "Profili Sottili") 
plot_button3 = Button(master = window, 
					command = met_trasf_conf, 
					height = 10, 
					width = 20, 
					text = "Jukowsky") 

# place the buttons 
# in main window 
plot_button.grid(row=1,column=1,sticky=N)
plot_button2.grid(row=1,column=2,sticky=N)
plot_button3.grid(row=1,column=3,sticky=N) 

#piazzare input di testo in cui mettere i dati
nacaa_var=StringVar()
nacaa_entry = ttk.Entry(window, width=7, textvariable=nacaa_var)
nacaa_entry.grid(column=1, row=2, sticky=(W, E))

U_inf_value=int()
curv_entry = ttk.Entry(window, width=7, textvariable=U_inf_value)
curv_entry.grid(column=2, row=2, sticky=(W, E))

AoA_value=float()
curv_entry = ttk.Entry(window, width=7, textvariable=AoA_value)
curv_entry.grid(column=3, row=2, sticky=(W, E))



nacaa_entry.focus()
# window.bind("<Return>", calculate())

#stile un po' più estetico
sv_ttk.use_dark_theme()
# run the gui 
window.mainloop()