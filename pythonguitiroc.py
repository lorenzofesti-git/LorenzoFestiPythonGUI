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


def plot1(): 
    #impostiamo i dati presi dalla GUI, e creaiamo un linspace
        # the figure that will contain the plot 
    fig = Figure(figsize = (5,4), dpi = 100)  
    # adding the subplot 
    plot1 = fig.add_subplot(111) 
    corda= 1
    x=np.array(linspace(0, float(corda), num=1000))
    naca=nacaa_var.get()
    #dividiamo naca cifra per cifra
    arr_naca=list(naca)
    naca_subtype = len(arr_naca)
    #normalizziamo per la corda
    xn=x/float(corda)
    xb=zeros(shape=(200,1))
    yb=zeros(shape=(200,1))
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


  
def calculate(*args):
    try:
        value = float(nacaa_var.get())
        AoA_value.set(float(value*0.1))
    except ValueError:
        pass

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