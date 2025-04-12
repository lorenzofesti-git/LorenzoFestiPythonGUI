
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


def nacanum(naca,x):
    corda=1
    arr_naca=list(naca)
    naca_subtype = len(arr_naca)
    xn=np.asarray(x/float(corda))
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
                deltY=max(yf)-min(yf)
                return xu,xl,yu,yl,xf,yf,ymediana,deltY
            
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
                deltY=max(yf)-min(yf)
                return xu,xl,yu,yl,xf,yf,ymediana,deltY
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
            deltY=max(yf)-min(yf)
            return xu,xl,yu,yl,xf,yf,ymediana,deltY
    else:
            return None  # Restituisce None se il numero non ha 4 o 5 arr_naca