#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

def parse_userparam() :
    from _UserParam import DataSet
    equation=DataSet["Probleme"]["equation"]
    
    T=DataSet["Domain"]["time"] 
    L = DataSet["Domain"]["length"]
    
    scheme=DataSet["Discretization"]["scheme"]
    F=DataSet["Discretization"]["fourier number"] 
    C=DataSet["Discretization"]["courant number"]         
    Nx=DataSet["Discretization"]["spatial resolution"] 
    version=DataSet["Discretization"]["version"]    
    
    animate=DataSet["Animation"]["active"]
    skip=DataSet["Animation"]["skip"] 

    set_init=DataSet["Init"]["active"]      
    init_loc=DataSet["Init"]["loc"]    
    init_shape=DataSet["Init"]["shape"] 
    init_sigma=DataSet["Init"]["sigma"] 

    set_source=DataSet["Source"]["active"]      
    source_loc=DataSet["Source"]["loc"]    
    source_type=DataSet["Source"]["type"] 
    source_sigma=DataSet["Source"]["sigma"] 
   
    a_0 = DataSet["Medium"]["diffusion_coeff"]
    c_0 = DataSet["Medium"]["celerity"]
    active_medium = DataSet["Medium"]["active_zone"]
    medium=DataSet["Medium"]["zone"] 
    slowness_factor=DataSet["Medium"]["slowness_factor"] 
    
    return equation, T, L, scheme, C, F, Nx, version, animate, skip, \
    set_init, init_loc, init_shape, init_sigma, \
    set_source, source_loc, source_type, source_sigma, \
    a_0, c_0, active_medium, medium, slowness_factor

def set_discretization(equation, L, Nx, a_0, c_0, C, F) :
    dx=(L/Nx)    
    
    if equation == "Diffusion":
        dt=F/a_0*dx**2
        
    elif equation in ("Advec","Wave"):
        dt=C*dx/c_0

    elif equation == "DiffAdvec":
        dt=np.min(C*dx/c_0,F/a_0*dx**2) 
        
    elif equation == "Burger":
        dt=F/a_0*dx**2

    elif equation == "inviscidBurger":
        dt=C*dx       
        
    else:
        raise ValueError('This equation name "%s" is wrong' % equation)     

    print('--> Grid discretization : ')
    print('- dx :',dx)
    print('- dt :',dt,"\n")
        
    return dx, dt
    
def set_casename(equation, T, L, scheme, C, F, dx, dt, \
    active_init, init_loc, init_shape, init_sigma, \
    active_source, source_loc, source_type, source_sigma, \
    a_0, c_0, active_medium, medium, slowness_factor):
    
    ADIM="";IMPULSION="";MEDIUM = ""
    
    if equation == 'Diffusion':
        ADIM = "a%s_F%s"%(a_0, F)
    elif equation == 'DiffAdvec':
        ADIM = "a%s_c%s_C%s_F%"%(a_0, c_0, C, F)
    else:
        ADIM = "c%s_C%s"%(c_0,C)
        
    if active_init == True :
        IMPULSION = "_init_%s_loc%s_sigma%s" % (init_shape, init_loc, init_sigma)
        
    if active_source == True :
        IMPULSION = "_source_%s_loc%s_sigma%s" % (source_shape, source_loc, source_sigma)
   
    if active_medium == True :
        MEDIUM = "-medium_%s_slowness_factor%s" % (medium, slowness_factor)

    casename = "%s_%s_T%s_L%s_dx%s_dt%s_%s%s%s" % (equation, scheme, T, L, dx, dt, ADIM, IMPULSION, MEDIUM)

    return casename

def set_window(equation) :
    if equation in ('Burger','inviscidBurger'):
        umin=-1.5; umax=1.5
    else:
        umin=-0.5; umax=1.5
    return umin, umax

def read_archive(equation, archive_name):

    array_names = np.load(archive_name)
    
    umin, umax = set_window(equation)
    n=0;  
    for array_name in array_names:
        if array_name != 't':
            if n == 0:
                plt.ion()
                lines = plt.plot(array_names['x'], array_names[array_name], 'r-')
                plt.axis([array_names['x'][0], array_names['x'][-1], umin, umax])
                plt.xlabel('x')
                plt.ylabel('u')
                plt.title("lecture de l'archive" + archive_name, fontsize=10)
                plt.legend(['t=%.3f' % array_names['t'][n-1]])
                plt.pause(0.0001) 
            else:
                # Update new solution
                lines[0].set_ydata(array_names[array_name])
                plt.legend(['t=%.3f' % array_names['t'][n-1]])
                plt.draw()
                plt.pause(0.0001)
            n=n+1