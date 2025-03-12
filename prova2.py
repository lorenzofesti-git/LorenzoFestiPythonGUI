from tkinter import *
from tkinter import ttk
import sv_ttk
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk) 
from matplotlib import *
from numpy import *
import numpy as np 

# list of squares

def plot1(): 
  
    # the figure that will contain the plot 
    fig = Figure(figsize = (5, 5), 
                 dpi = 100) 
  
    # adding the subplot 
    plot1 = fig.add_subplot(111) 
    plot1.set_xlim(0,100)
    #RANGE PARTE SEMPRE DA ZERO! PYTHON INDICIZZA A PARTIRE DA ZERO!!

    corda= chord_value.get()
    x=linspace(0, int(corda), num=1000)
    naca=nacaa_var.get()
    arr_naca=list(naca)
    naca_subtype = len(arr_naca)
    if naca_subtype == 4:
        m=(int(arr_naca[0]))/100
        p=int(arr_naca[1])/10
        t=(10*int(arr_naca[2])+int(arr_naca[3]))/100
        if m==0 and p==0:
            yt=(t/0.2)((0.2969*x(0.5))-(0.126*x)-(0.3516*x2)+(0.2843*x3)-(0.1015*x*4))
            xu=x 
            xl=x
            yu=yt
            yl=-yt
            
        else:
            yt=(t/0.2)((0.2969*x(0.5))-(0.126*x)-(0.3516*x2)+(0.2843*x3)-(0.1015*x*4))
    elif naca_subtype == 5:
        m=1.5*int(arr_naca[0])
        p=0.5*(10*int(arr_naca[1])+int(arr_naca[2]))
        t=10*int(arr_naca[3])+int(arr_naca[4])
    else:
        return None  # Restituisce None se il numero non ha 4 o 5 arr_naca
    plot1.plot(xu,yu,'b-',xl,yl,'b-') 
  
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

# button that displays the plot 
plot_button = Button(master = window, 
					command = plot1, 
					height = 10, 
					width = 20, 
					text = "Jukowsky") 
plot_button2 = Button(master = window, 
					command = plot1, 
					height = 10, 
					width = 20, 
					text = "Vorticità puntiformi") 
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


nacaa_var=StringVar()
nacaa_entry = ttk.Entry(window, width=7, textvariable=nacaa_var)
nacaa_entry.grid(column=2, row=2, sticky=(W, E))




chord_value=StringVar()
curv_entry = ttk.Entry(window, width=7, textvariable=chord_value)
curv_entry.grid(column=3, row=2, sticky=(W, E))

nacaa_entry.focus()
window.bind("<Return>", calculate())

# w=Spinbox()  per modificare i valori

#stile un po' più estetico
sv_ttk.use_dark_theme()
# run the gui 
window.mainloop()