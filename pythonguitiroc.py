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
from import1 import *
from import5 import *
from import6 import *
from nacanumber import *

def plot1(): 
    #impostiamo i dati presi dalla GUI, e creaiamo un linspace
        # the figure that will contain the plot 
    fig = Figure(figsize = (5,4), dpi = 100)  
    # adding the subplot 
    plot1 = fig.add_subplot(111) 
    corda= 1
    x=np.array(linspace(0, float(corda), num=1000))
    naca=nacaa_var.get()
    xb=zeros(shape=(200,1))
    yb=zeros(shape=(200,1))
    nacanum(naca,x)
    #  IMPORTA QUI DENTRO LA FUNZIONE PER CREARE PANNELLI


    # definiamo spazio del plot

        # creating the Tkinter canvas 
    # containing the Matplotlib figure 
    canvas = FigureCanvasTkAgg(fig, 
                                master = window)   
    canvas.get_tk_widget().grid(row=5,column=1,columnspan=3,rowspan=20,sticky=N)
        # here: plot suff to your fig

    canvas.draw()
    
    plot1.set_xlim(-0.5,1.2)
    plot1.set_ylim(-max(xu)*1.2,max(xu)*1.2)
    plot1.plot(xu,yu,'k-',xl,yl,'k-',x,ymediana,'r-',xb,yb,'2') 

#
# def calculate(*args):
#     try:
#         value = float(nacaa_var.get())
#         AoA_value.set(float(value*0.1))
#     except ValueError:
#         pass
#
# the main Tkinter window 
window = Tk() 
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
					command = plot2, 
					height = 10, 
					width = 20, 
					text = "Profili Sottili") 
plot_button3 = Button(master = window, 
					command = plot3, 
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
nacaa_entry.grid(column=1, row=2, sticky=(W, E))

AoA_value=StringVar()
curv_entry = ttk.Entry(window, width=7, textvariable=AoA_value)
curv_entry.grid(column=3, row=2, sticky=(W, E))

U_inf_value=StringVar()
curv_entry = ttk.Entry(window, width=7, textvariable=U_inf_value)
curv_entry.grid(column=2, row=2, sticky=(W, E))

nacaa_entry.focus()
window.bind("<Return>", calculate())

#stile un po' più estetico
sv_ttk.use_dark_theme()
# run the gui 
window.mainloop()