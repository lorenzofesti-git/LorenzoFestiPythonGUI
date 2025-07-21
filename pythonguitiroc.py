from tkinter import *
import sv_ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.widgets import Slider
from matplotlib import *
from plotly.figure_factory import create_streamline
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
from lasthope import *


def met_panels():
    # impostiamo i dati presi dalla GUI, e creaiamo un linspace
    # the figure that will contain the plot
    fig = Figure(figsize=(5, 4), dpi=100)
    # adding the subplot
    met_panels = fig.add_subplot(111)
    corda = 1
    x = np.array(linspace(0, float(corda), num=100))
    naca = nacaa_entry.get()
    U_inf_value = curv_entry1.get()
    AoA_value = curv_entry2.get()

    xu, xl, yu, yl, xf, yf, ymediana, ytrack, dydx = nacanum(naca, x)

    xb, yb, VxC, VyC, XX, YY = pannelli(
        float(U_inf_value), float(AoA_value), xf, yf, xu, yu, xl, yl
    )
    type = 1

    # definiamo spazio del plot

    # creating the Tkinter canvas
    # containing the Matplotlib figure
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.get_tk_widget().grid(row=5, column=1, columnspan=3, rowspan=20, sticky=N)
    # here: plot suff to your fig

    canvas.draw()

    met_panels.set_xlim(-0.2, 1.2)
    met_panels.set_ylim(-0.5, 0.5)
    met_panels.plot(xu, yu, "k-", xl, yl, "k-")
    met_panels.streamplot(XX, YY, VxC, VyC)
    met_panels.grid()
    # ,xb,yb,'2'


def met_prof_sott():
    # create the figure that will contain the plot
    fig = Figure(figsize=(5, 4), dpi=100)
    # adding the subplot
    met_panels = fig.add_subplot(111)
    corda = 1
    x = np.array(linspace(0, float(corda), num=100))
    naca = nacaa_var.get()
    U_inf_value = curv_entry1.get()
    AoA_value = curv_entry2.get()
    xu, xl, yu, yl, xf, yf, ymediana, yt, dydx = nacanum(naca, x)
    X, Y, U, V, x_colloc, gamma, cl, cm_c4, cp = calculate_thin_airfoil(
        x, ymediana, int(U_inf_value), int(AoA_value), dydx, xu, yu, xl, yl
    )

    type = 2
    # gam,cl,cm_c4,x_cp,u_x,v_x,XX, YY = solve_thin_airf_theory(x,ymediana,U_inf_value,AoA_value)

    # define space to plot into

    # creating the Tkinter canvas
    # containing the Matplotlib figure
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.get_tk_widget().grid(row=5, column=1, columnspan=3, rowspan=20, sticky=N)
    # here: plot suff to your fig

    canvas.draw()
    met_panels.set_xlim(-0.2, 1.2)
    met_panels.set_ylim(-1.5, 1.5)
    # met_panels.plot(xu,yu,'g-',xl,yl,'r-')
    met_panels.plot(x, ymediana, "k-")
    met_panels.streamplot(X, Y, U, V)
    met_panels.grid()


def met_trasf_conf():
    # the figure that will contain the plot
    fig = Figure(figsize=(5, 4), dpi=100)
    # adding the subplot

    met_panels = fig.add_subplot(111)
    corda = 1
    x = np.array(linspace(0, float(corda), num=100))
    # naca=nacaa_var.get()
    U_inf_value = curv_entry1.get()
    AoA_value = curv_entry2.get()
    k = nacaa_var.get()
    c = int(k) / 100

    f, J, zcirc, zair = vectorized_version(U_inf_value, AoA_value, c)
    # creating the Tkinter canvas
    # containing the Matplotlib figure
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.get_tk_widget().grid(row=5, column=1, columnspan=3, rowspan=20, sticky=N)
    # here: plot suff to your fig
    met_panels.contour(np.real(J), np.imag(J), np.imag(f), 100)
    met_panels.fill(np.real(zair), np.imag(zair), "r")
    met_panels.set_aspect("equal")
    met_panels.set_xlim([-2.5, 2.5])
    met_panels.set_ylim([-2.5, 2.5])
    met_panels.set_title("Joukowski Airfoil")
    met_panels.grid(True, alpha=0.3)

    plt.show()


# the main Tkinter window
window = Tk()
window.columnconfigure(3)
window.rowconfigure(3)

# setting the title
window.title("Plotting")

# dimensions of the main window

# buttons that displays the 3 plots:

# panel method
plot_button = Button(
    master=window, command=met_panels, height=10, width=20, text="Metodo a Pannelli"
)

# thion airfoil theory
plot_button2 = Button(
    master=window, command=met_prof_sott, height=10, width=20, text="Profili Sottili"
)
# jouw. transf.
plot_button3 = Button(
    master=window, command=met_trasf_conf, height=10, width=20, text="Jukowsky"
)

# add the save function
plot_button4 = Button(
    master=window, command=met_trasf_conf, height=10, width=20, text="Save results"
)

# place the buttons
# in main window
plot_button.grid(row=1, column=1, sticky=N)
plot_button2.grid(row=1, column=2, sticky=N)
plot_button3.grid(row=1, column=3, sticky=N)
plot_button4.grid(column=1, row=10, columnspan=20, sticky=(S))

# piazzare input di testo in cui mettere i dati
text = Label(window, text="numero NACA")
text.grid(column=0, row=2, sticky=(W, E))

nacaa_var = StringVar()
nacaa_entry = ttk.Entry(window, width=7, textvariable=nacaa_var)
nacaa_entry.grid(column=1, row=2, sticky=(W, E))

text = Label(window, text="modulo velocità")
text.grid(column=2, row=2, sticky=(W, E))

U_inf_value = int()
curv_entry1 = ttk.Entry(window, width=7, textvariable=U_inf_value)
curv_entry1.grid(column=3, row=2, sticky=(W, E))

text = Label(window, text="modulo angolo di attacco")
text.grid(column=4, row=2, sticky=(W, E))

AoA_value = float()
curv_entry2 = ttk.Entry(window, width=7, textvariable=AoA_value)
curv_entry2.grid(column=5, row=2, sticky=(W, E))


nacaa_entry.focus()
# window.bind("<Return>", calculate())

# stile un po' più estetico
sv_ttk.use_dark_theme()
# run the gui
window.mainloop()
