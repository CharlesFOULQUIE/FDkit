#!/usr/bin/env python

import time, os
import numpy as np
from Modules import NumericalScheme

def make_hash(I,V,f,c,U_0,U_L,L,dx,dt,C,F,T,equation,scheme):
    import hashlib, inspect
    data = inspect.getsource(I) + '_' + inspect.getsource(V) + \
           '_' + inspect.getsource(f) + '_' + str(c) + '_' + \
           ('None' if U_0 is None else inspect.getsource(U_0)) + \
           ('None' if U_L is None else inspect.getsource(U_L)) + \
           '_' + str(L) + '_' + str(dx) + '_' + str(dt) + '_' + \
           str(C) + '_' + str(T) + '_' + str(equation) + '_' + str(scheme)
    hashed_input = hashlib.sha1(data.encode('utf-8')).hexdigest()
    return hashed_input



def Initialize(f, I, V, U_0, U_L, c, x, Nx, t, version):
    if f is None or f == 0:
        f = (lambda x, t: 0) if version == 'scalar' else \
            lambda x, t: np.zeros(x.shape)
    if I is None or I == 0:
        I = (lambda x: 0) if version == 'scalar' else \
            lambda x: np.zeros(x.shape)
    if V is None or V == 0:
        V = (lambda x: 0) if version == 'scalar' else \
            lambda x: np.zeros(x.shape)
    if U_0 is not None:
        if isinstance(U_0, (float,int)) and U_0 == 0:
            U_0 = lambda t: 0
    if U_L is not None:
        if isinstance(U_L, (float,int)) and U_L == 0:
            U_L = lambda t: 0
            
    # Make c(x) available as array
    if isinstance(c, (float,int)):
        c = np.zeros(x.shape) + c
    elif callable(c):
        # Call c(x) and fill array c
        c_ = np.zeros(x.shape)
        for i in range(Nx+1):
            c_[i] = c(x[i])
        c = c_
        
    return f, I, V, U_0, U_L, c
    
def TimeLoop(equation, scheme, u, u_n, u_nm1, x, t, dx, dt, c, U_0, U_L, Nx, Nt, V, f, C, F, version, user_action, bc_type):
    It = range(0, Nt+1)
    
    if equation == "Wave":
    
        if scheme == "CD":
            for n in It[0:-1]:
                u, u_n, u_nm1 = NumericalScheme.wave1D(u, u_n, u_nm1, n, x, t, dx, dt, c, U_0, U_L, Nx, V, f, version, bc_type) # --- Compute
                if user_action is not None: # --- Plot solution
                    if user_action(u, x, t, n+1):
                        break
                u_nm1, u_n, u = u_n, u, u_nm1 # --- Update data structures for next step
        else:
            raise ValueError('{} not implemented'.format(scheme)) 
                
    elif equation == "Diffusion":
        if scheme == "FE_CD":
            for n in It[0:-1]:
                u, u_n = NumericalScheme.diffusion1D_FE_CD(u, u_n, Nx, F, f, n, dt, x, t) # --- Compute
                if user_action is not None: # --- Plot solution
                    if user_action(u, x, t, n+1):
                        break
                u, u_n = u_n, u # --- Update data structures for next step
        else:
            raise ValueError('{} not implemented'.format(scheme)) 

    elif equation == "DiffAdvec":
        if scheme == "FE_LW":
            for n in It[0:-1]:
                u, u_n = NumericalScheme.Diff_advec1D_FE_LW(u, u_n, Nx, C, F, f, n, dt, x, t) # --- Compute
                if user_action is not None: # --- Plot solution
                    if user_action(u, x, t, n+1):
                        break
                u, u_n = u_n, u # --- Update data structures for next step 
        else:
            raise ValueError('{} not implemented'.format(scheme))   

    elif equation == "Burger":
        if scheme == "FE_LW":
            for n in It[0:-1]:
                u, u_n = NumericalScheme.Burger_FE_LW(u, u_n, Nx, C, F, f, n, dt, dx, x, t) # --- Compute
                if user_action is not None: # --- Plot solution
                    if user_action(u, x, t, n+1):
                        break
                u, u_n = u_n, u # --- Update data structures for next step 
        else:
            raise ValueError('{} not implemented'.format(scheme))   

    elif equation == "inviscidBurger":
        if scheme == "FE_LW":
            for n in It[0:-1]:
                u, u_n = NumericalScheme.inviscidBurger_FE_LW(u, u_n, Nx, C, f, n, dt, dx, x, t) # --- Compute
                if user_action is not None: # --- Plot solution
                    if user_action(u, x, t, n+1):
                        break
                u, u_n = u_n, u # --- Update data structures for next step 
        else:
            raise ValueError('{} not implemented'.format(scheme)) 
 
    elif equation == "Advec": 
    
        if scheme == "Leapfrog": 
            for n in It[0:-1]:
                u, u_n, u_nm1 = NumericalScheme.advec1D_Leapfrog(u, u_n, u_nm1, n, Nx, C, f, dt, x, t, bc_type) # --- Compute
                if user_action is not None: # --- Plot solution
                    if user_action(u, x, t, n+1):
                        break
                u_nm1, u_n, u = u_n, u, u_nm1 # --- Update data structures for next step
                

        elif scheme == "LW":
            for n in It[0:-1]:
                u, u_n = NumericalScheme.advec1D_LW(u, u_n, Nx, C, f, n, dt, x, t, bc_type) # --- Compute
                if user_action is not None: # --- Plot solution
                    if user_action(u, x, t, n+1):
                        break
                u, u_n = u_n, u # --- Update data structures for next step

        elif scheme == "LF":
            for n in It[0:-1]:
                u, u_n = NumericalScheme.advec1D_LF(u, u_n, Nx, C, f, n, dt, x, t, bc_type) # --- Compute
                if user_action is not None: # --- Plot solution
                    if user_action(u, x, t, n+1):
                        break
                u, u_n = u_n, u # --- Update data structures for next step

        elif scheme == "UP":
            for n in It[0:-1]:
                u, u_n = NumericalScheme.advec1D_UP(u, u_n, Nx, C, f, n, dt, x, t, bc_type) # --- Compute
                if user_action is not None: # --- Plot solution
                    if user_action(u, x, t, n+1):
                        break
                u, u_n = u_n, u # --- Update data structures for next step

        elif scheme == "FE_CD": # Forward Euler scheme in time and centered differences in space (FECS).
            for n in It[0:-1]:
                u, u_n = NumericalScheme.advec1D_FE_CD(u, u_n, Nx, C, f, n, dt, x, t, bc_type) # --- Compute
                if user_action is not None: # --- Plot solution
                    if user_action(u, x, t, n+1):
                        break
                u, u_n = u_n, u # --- Update data structures for next step
        else:
            raise ValueError('{} not implemented'.format(scheme))   

    else:
        raise ValueError('{} not implemented'.format(equation))        

    return 0

def solver(
    equation, scheme, 
    I, V, f, c, U_0, U_L, L, dt, dx, C, F, T,
    user_action=None, version='scalar', bc_type=None):

    # --- Compute time and space grid ---
    Nt = int(round(T/dt))
    t = np.linspace(0, Nt*dt, Nt+1)      # Mesh points in time
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)          # Mesh points in space
    
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    # --- Initialize input data f, I, V, U_0, U_L if None or 0 ---
    f, I, V, U_0, U_L, c = Initialize(f, I, V, U_0, U_L, c, x, Nx, t, version)

    # --- Make hash of all input data ---
    hashed_input=make_hash(I,V,f,c,U_0,U_L,L,dx,dt,C,F,T,equation,scheme)
    
    # Check if simulation is already run
    check = True
    if check == True :
        if os.path.isfile(os.path.join(user_action.get_res_directory(),'.' + hashed_input + '_archive.npz')):
            return -1, hashed_input

    # --- Allocate memomry for solutions ---
    u     = np.zeros(Nx+1)   # Solution array at new time level
    u_n   = np.zeros(Nx+1)   # Solution at 1 time level back
    u_nm1 = np.zeros(Nx+1)   # Solution at 2 time levels back    
    
    # --- Load initial condition into u_n ---
    for i in range(0,Nx+1):
        u_n[i] = I(x[i])
            
    # --- Plot and Store Solution ---
    if user_action is not None:
        user_action(u_n, x, t, 0)

    print("\n")
    print("----- Begin Time loop -----")
    print("\n")

    # --- Time loop ---
    t0 = time.time()  # CPU time measurement
    TimeLoop(equation, scheme, u, u_n, u_nm1, x, t, dx, dt, c, U_0, U_L, Nx, Nt, V, f, C, F, version, user_action, bc_type)
    cpu_time = time.time() - t0
    
    print("\n")
    print("----- End Time loop -----")
    print("\n")
    
    return cpu_time, hashed_input


