#!/usr/bin/env python

def diffusion1D_FE_CD(u, u_n, Nx, F, f, n, dt, x, t): 
    ''' Forward Euler scheme in time '''
    for i in range(1, Nx):
        u[i] = u_n[i] + F*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + dt*f(x[i], t[n])

    # Insert boundary conditions
    u[0] = 0;  u[Nx] = 0;
    
    return u, u_n
    
def Diff_advec1D_FE_LW(u, u_n, Nx, C, F, f, n, dt, x, t):
    ''' The Lax-Wendroff method '''
    
    return u, u_n
    

def advec1D_Leapfrog(u, u_n, u_nm1, n, Nx, C, f, dt, x, t):
    ''' Leapfrog in time, centered differences in space '''
    
    return u, u_n, u_nm1


def advec1D_LW(u, u_n, Nx, C, f, n, dt, x, t):
    ''' The Lax-Wendroff method '''

    
    return u, u_n
    
    
def advec1D_LF(u, u_n, Nx, C, f, n, dt, x, t):
    ''' The Lax–Friedrichs method '''

    
    return u, u_n
    

def advec1D_UP(u, u_n, Nx, C, f, n, dt, x, t):
    ''' Forward Euler scheme in time and UPwind differences in space '''
    periodic_bc=True

    
    return u, u_n


def advec1D_FE_CD(u, u_n, Nx, C, f, n, dt, x, t): 
    ''' Forward Euler scheme in time and centered differences in space (FECS)  '''
    periodic_bc=True

    if periodic_bc:
        i = 0
        u[i] = u_n[i] - 0.5*C*(u_n[i+1] - u_n[Nx]) + dt*f(x[i], t[n])
        u[Nx] = u[0]
        
    for i in range(1, Nx):
        u[i] = u_n[i] - 0.5*C*(u_n[i+1] - u_n[i-1]) + dt*f(x[i], t[n])

    if not periodic_bc:
        u[0] = 0
    
    return u, u_n



def wave1D(u, u_n, u_nm1, n, x, t, dx, dt, c, U_0, U_L, Nx, V, f, version):

    Ix = range(0, Nx+1)
    C2 = (dt/dx)**2; 
    dt2 = dt*dt
    
    if n == 0:
        # --- Special formula for the first step ---
        for i in Ix[1:-1]:
            u[i] = u_n[i] + dt*V(x[i]) + \
            0.5*C2*(0.5*(c[i]**2 + c[i+1]**2)*(u_n[i+1] - u_n[i]) - \
                    0.5*(c[i]**2 + c[i-1]**2)*(u_n[i] - u_n[i-1])) + \
            0.5*dt2*f(x[i], t[0])
    
        i = Ix[0]
        if U_0 is None:
            # Set boundary values (x=0: i-1 -> i+1 since u[i-1]=u[i+1]
            # when du/dn = 0, on x=L: i+1 -> i-1 since u[i+1]=u[i-1])
            ip1 = i+1
            im1 = ip1  # i-1 -> i+1
            u[i] = u_n[i] + dt*V(x[i]) + \
                0.5*C2*(0.5*(c[i]**2 + c[ip1]**2)*(u_n[ip1] - u_n[i])  - \
                        0.5*(c[i]**2 + c[im1]**2)*(u_n[i] - u_n[im1])) + \
            0.5*dt2*f(x[i], t[0])
        else:
            u[i] = U_0(dt)
    
        i = Ix[-1]
        if U_L is None:
            im1 = i-1
            ip1 = im1  # i+1 -> i-1
            u[i] = u_n[i] + dt*V(x[i]) + \
                0.5*C2*(0.5*(c[i]**2 + c[ip1]**2)*(u_n[ip1] - u_n[i])  - \
                        0.5*(c[i]**2 + c[im1]**2)*(u_n[i] - u_n[im1])) + \
            0.5*dt2*f(x[i], t[0])
        else:
            u[i] = U_L(dt)
    else:
        # Update all inner points
        if version == 'scalar':
            for i in Ix[1:-1]:
                u[i] = - u_nm1[i] + 2*u_n[i] + \
                    C2*(0.5*(c[i]**2 + c[i+1]**2)*(u_n[i+1] - u_n[i])  - \
                        0.5*(c[i]**2 + c[i-1]**2)*(u_n[i] - u_n[i-1])) + \
                dt2*f(x[i], t[n])
    
        elif version == 'vectorized':
            u[1:-1] = - u_nm1[1:-1] + 2*u_n[1:-1] + \
            C2*(0.5*(c[1:-1]**2 + c[2:]**2)*(u_n[2:] - u_n[1:-1]) -
                0.5*(c[1:-1]**2 + c[:-2]**2)*(u_n[1:-1] - u_n[:-2])) + \
            dt2*f(x[1:-1], t[n])
        else:
            raise ValueError('version=%s' % version)
    
        # Insert boundary conditions
        i = Ix[0]
        if U_0 is None:
            # Set boundary values
            # x=0: i-1 -> i+1 since u[i-1]=u[i+1] when du/dn=0
            # x=L: i+1 -> i-1 since u[i+1]=u[i-1] when du/dn=0
            ip1 = i+1
            im1 = ip1
            u[i] = - u_nm1[i] + 2*u_n[i] + \
                C2*(0.5*(c[i]**2 + c[ip1]**2)*(u_n[ip1] - u_n[i])  - \
                    0.5*(c[i]**2 + c[im1]**2)*(u_n[i] - u_n[im1])) + \
            dt2*f(x[i], t[n])
        else:
            u[i] = U_0(t[n+1])
    
        i = Ix[-1]
        if U_L is None:
            im1 = i-1
            ip1 = im1
            u[i] = - u_nm1[i] + 2*u_n[i] + \
                C2*(0.5*(c[i]**2 + c[ip1]**2)*(u_n[ip1] - u_n[i])  - \
                    0.5*(c[i]**2 + c[im1]**2)*(u_n[i] - u_n[im1])) + \
            dt2*f(x[i], t[n])
        else:
            u[i] = U_L(t[n+1])
        
    return u, u_n, u_nm1