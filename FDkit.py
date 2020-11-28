#!/usr/bin/env python

''' Auteur : Charles Foulquié & Hans Petter Langtangen
    License : CC BY Attribution 4.0 license 
    Historique : Modification du projet GitHub hplgit écrit par Hans Petter Langtangen '''

import time, glob, shutil, os
import numpy as np
from Modules import PlotAndStoreSolution, Solver, Utils

def set_celerity(c_0, medium, slowness_factor, active):
    global c
    if active == True :
        def c(x):
            return c_0/slowness_factor \
                if medium[0] <= x <= medium[1] else c_0
    else :
        def c(x):
            return c_0

def set_source(active, stype,x0,sigma):
    global f
    if active == True :
        if stype in ('harmonic_gaussian_peak'):
            def f(x,t):
                return 100*np.exp(-0.5*((x-x0)/sigma)**2)*np.cos(2*np.pi*2.0*t) 
        elif stype in ('progressive_gaussian_peak'):
            def f(x,t):
                return np.exp(-0.5*((x-x0)/sigma)**2) \
                    if t <= 0.25 else 0
        else:
            raise ValueError('Wrong source setting="%s"' % stype)
    else :
        def f(x,t):
            return 0
    return 0
    
def set_bc(bc_type):
    global U_0
    global U_L
    f=10

    if bc_type['left'] in ('harmonic'):
        def U_0(t):
            return np.cos(2*np.pi*f*t) 
    else :
        def U_0(t):
            return 0    
    
    if bc_type['right'] in ('harmonic'):
        def U_L(t):
            return np.cos(2*np.pi*f*t) 
    else :
        def U_L(t):
            return 0
    return 0

    
def set_init(active, ishape, x0, sigma):
    global I

    if active == True :
        if ishape in ('gaussian','Gaussian'):
            def I(x):
                return np.exp(-0.5*((x-x0)/sigma)**2)
        elif ishape == 'plug':
            def I(x):
                return np.where(abs(x-x0) > sigma, 0, 1)
                #return 0 if abs(x-x0) > sigma else 1
        elif ishape == 'cosinehat':
            def I(x):
                # One period of a cosine
                w = 2
                a = w*sigma
                return np.where((x0 - a <= x) & (x <= x0 + a), 0.5*(1 + np.cos(np.pi*(x-x0)/a)), 0)
                #return 0.5*(1 + np.cos(np.pi*(x-x0)/a)) \
                #    if x0 - a <= x <= x0 + a else 0
    
        elif ishape == 'half-cosinehat':
            def I(x):
                # Half a period of a cosine
                w = 4
                a = w*sigma
                return np.where((x0 - 0.5*a <= x) & (x <= x0 + 0.5*a), np.cos(np.pi*(x-x0)/a), 0)
                #return np.cos(np.pi*(x-init_loc)/a) \
                #    if init_loc - 0.5*a <= x <= init_loc + 0.5*a else 0
        
        elif ishape in ('sinus'):
            def I(x):
                return np.sin(2*np.pi*x)
        else:
            raise ValueError('Wrong peak shape="%s"' % init_shape)
    else :
        def I(x):
            return 0
    return 0
    


def main():

    equation, T, L, scheme, C, F, Nx, version, animate, skip, \
    active_init, init_loc, init_shape, init_sigma, \
    active_source, source_loc, source_type, source_sigma, \
    a_0, c_0, active_medium, medium, \
    slowness_factor, bc_type= Utils.parse_userparam()
    
    # set grid
    dx, dt = Utils.set_discretization(equation, L, Nx, a_0, c_0, C, F)
    
    # set casename
    casename=Utils.set_casename(equation, T, L, scheme, C, F, dx, dt, \
    active_init, init_loc, init_shape, init_sigma, \
    active_source, source_loc, source_type, source_sigma, \
    a_0, c_0, active_medium, medium, slowness_factor)
  
    print("   _____   .___ __    .__   __    ") 
    print(" _/ ____\__| _/|  | __|__|_/  |_  ")
    print(" \   __\/ __ | |  |/ /|  |\   __\ ")
    print("  |  | / /_/ | |    < |  | |  |   ")
    print("  |__| \____ | |__|_ \|__| |__|   ")
    print("            \/      \/            ")
    print("------** Problem to solve **------")
    print("Equation : ", equation,"\n")
    
    print("------** Domain **------")
    print("Time : ", T)
    print("Length : ", L, "\n")
    
    print("------** Discretization **------")
    print("Scheme : ", scheme)
    print("Fourier number : ", F)
    print("Courant number : ", C)
    print("Spatial resolution : ", L)
    print("version : ", version,"\n")
    
    print("------** Init **------")
    print("set init : ", active_init)
    print("init peak loc : ", init_loc)
    print("init peak shape : ", init_shape)
    print("init peak width : ", init_sigma)
    
    print("------** Source **------")
    print("set source : ", active_source)
    print("source peak loc : ", source_loc)
    print("source peak shape : ", source_type)
    print("source peak width : ", source_sigma,"\n")
    
    print("------** Medium **------")
    print("diffusion coeff : ", a_0)
    print("celerity : ", c_0)
    print("active_zone : ", active_medium)
    print("zone : ", medium)
    print("slowness_factor : ", slowness_factor,"\n")


    # Set test case configuration
    set_init(active_init, init_shape, init_loc, init_sigma)
    set_source(active_source, source_type, source_loc, source_sigma) 
    set_celerity(c_0, medium, slowness_factor, active_medium)
    set_bc(bc_type)
    
    umin, umax = Utils.set_window(equation, bc_type)
    
    # Set animation style
    if active_medium == True :
        action = PlotAndStoreSolution.PlotMediumAndSolution(
            medium, casename=casename, umin=umin, umax=umax,
            skip=skip, screen_movie=animate, filename='tmpdata', init=I)
    else :
        action = PlotAndStoreSolution.PlotAndStoreSolution(
            casename=casename, umin=umin, umax=umax,
            skip=skip, screen_movie=animate, filename='tmpdata', init=I)
    
    # Solve problem
    cpu, hashed_input = Solver.solver(
        equation, scheme, I=I, V=None, f=f, 
        c=c, U_0=U_0, U_L=U_L, L=L, dt=dt, 
        dx=dx, C=C, F=F, T=T, user_action=action, 
        version=version, bc_type=bc_type)

    if cpu > 0:  # did we generate new data?
        print('hashed_input :', hashed_input)
        action.close_file(hashed_input)
        action.make_movie_file()
        print('*** cpu time :', cpu)
    else:
        print('Des données ont déjà été généréess avec ce jeu de paramètre !!')
        print("Souhaitez vous lire l'archive ?")
        user_choice = input()
        if user_choice in ('oui','1','OUI','Oui','ouii'):
            Utils.read_archive(equation, os.path.join(action.get_res_directory(),'.' + hashed_input + '_archive.npz'))
        else:
            print("ok tant pi")

if __name__ == '__main__':
    main();