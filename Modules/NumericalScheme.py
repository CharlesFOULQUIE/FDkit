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
    
def advec1D_Leapfrog(u, u_n, u_nm1, n, Nx, C, f, dt, x, t, bc_type):
    ''' Leapfrog in time, centered differences in space '''
    periodic_bc = False
    if bc_type['left'] == bc_type['right'] and bc_type['right'] == 'periodic':
        periodic_bc = True

    if n == 0:
        # Use upwind for first step
        if periodic_bc:
            i = 0
            #u[i] = u_n[i] - C*(u_n[i] - u_n[Nx-1])
            u_n[i] = u_n[Nx]
        for i in range(1, Nx+1):
            u[i] = u_n[i] - C*(u_n[i] - u_n[i-1]) + dt*f(x[i], t[0])
    else:
        if periodic_bc:
            i = 0
            # Must have this,
            u[i] = u_nm1[i] - C*(u_n[i+1] - u_n[Nx-1]) + dt*f(x[i], t[n])
            # not this:
            #u_n[i] = u_n[Nx]
        for i in range(1, Nx):
            u[i] = u_nm1[i] - C*(u_n[i+1] - u_n[i-1]) + dt*f(x[i], t[n])
        if periodic_bc:
            u[Nx] = u[0]
    if not periodic_bc:
        u[0] = 0
        u[Nx] = 0
    
    return u, u_n, u_nm1


def advec1D_LW(u, u_n, Nx, C, f, n, dt, x, t, bc_type):
    ''' The Lax-Wendroff method '''
    periodic_bc = False
    if bc_type['left'] == bc_type['right'] and bc_type['right'] == 'periodic':
        periodic_bc = True

    if periodic_bc:
        i = 0
        # Must have this,
        u[i] = u_n[i] - 0.5*C*(u_n[i+1] - u_n[Nx-1]) + \
               0.5*C**2*(u_n[i+1] - 2*u_n[i] + u_n[Nx-1]) + dt*f(x[i], t[n])
        # not this:
        #u_n[i] = u_n[Nx]
    for i in range(1, Nx):
        u[i] = u_n[i] - 0.5*C*(u_n[i+1] - u_n[i-1]) + \
               0.5*C**2*(u_n[i+1] - 2*u_n[i] + u_n[i-1]) + dt*f(x[i], t[n])

    if periodic_bc:
        u[Nx] = u[0]
        
    if not periodic_bc:
        u[0] = 0
    
    return u, u_n
    
    
def advec1D_LF(u, u_n, Nx, C, f, n, dt, x, t, bc_type):
    ''' The Lax–Friedrichs method '''
    periodic_bc = False
    if bc_type['left'] == bc_type['right'] and bc_type['right'] == 'periodic':
        periodic_bc = True

    if periodic_bc:
        i = 0
        # Must have this,
        u[i] = 0.5*(u_n[i+1] + u_n[Nx-1]) - 0.5*C*(u_n[i+1] - u_n[Nx-1]) + dt*f(x[i], t[n])
        # not this:
        #u_n[i] = u_n[Nx]
    for i in range(1, Nx):
        u[i] = 0.5*(u_n[i+1] + u_n[i-1]) - 0.5*C*(u_n[i+1] - u_n[i-1])  + dt*f(x[i], t[n])

    if periodic_bc:
        u[Nx] = u[0]
        
    if not periodic_bc:
        u[0] = 0
    
    return u, u_n
    

def advec1D_UP(u, u_n, Nx, C, f, n, dt, x, t):
    ''' Forward Euler scheme in time and UPwind differences in space '''
    periodic_bc=True

    
    return u, u_n


def advec1D_FE_CD(u, u_n, Nx, C, f, n, dt, x, t, bc_type): 
    ''' Forward Euler scheme in time and centered differences in space (FECS)  '''
    periodic_bc = False
    if bc_type['left'] == bc_type['right'] and bc_type['right'] == 'periodic':
        periodic_bc = True

    if periodic_bc == True:
        i = 0
        u[i] = u_n[i] - 0.5*C*(u_n[i+1] - u_n[Nx]) + dt*f(x[i], t[n])
        u[Nx] = u[0]
        
    for i in range(1, Nx):
        u[i] = u_n[i] - 0.5*C*(u_n[i+1] - u_n[i-1]) + dt*f(x[i], t[n])

    if not periodic_bc:
        u[0] = 0
    
    return u, u_n



def wave1D(u, u_n, u_nm1, n, x, t, dx, dt, c, U_0, U_L, Nx, V, f, version, bc_type):

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
        if bc_type['left'] == 'reflective':
            # Set boundary values
            # x=0: i-1 -> i+1 since u[i-1]=u[i+1] when du/dn=0
            # x=L: i+1 -> i-1 since u[i+1]=u[i-1] when du/dn=0
            ip1 = i+1
            im1 = ip1
            u[i] = - u_nm1[i] + 2*u_n[i] + \
                C2*(0.5*(c[i]**2 + c[ip1]**2)*(u_n[ip1] - u_n[i])  - \
                    0.5*(c[i]**2 + c[im1]**2)*(u_n[i] - u_n[im1])) + \
            dt2*f(x[i], t[n])
            
        elif bc_type['left'] == 'non-reflective':
            # Set boundary values
            # x=0: i-1 -> i+1 since u[i-1]=u[i+1] when du/dn=0
            # x=L: i+1 -> i-1 since u[i+1]=u[i-1] when du/dn=0
            ip1 = i+1
            im1 = ip1
            u[i] = - u_nm1[i] + 2*u_n[i] + \
                C2*(0.5*(c[i]**2 + c[ip1]**2)*(u_n[ip1] - u_n[i])  - \
                    0.5*(c[i]**2 + c[i]**2)*(u_n[i] - u_nm1[i])) + \
            dt2*f(x[i], t[n])
            
        elif bc_type['left'] == 'harmonic':
            u[i] = U_0(t[n+1])
            
        else :
            raise ValueError('bc left error : {} not implemented'.format(bc_type['left'])) 
    
        i = Ix[-1]
        if bc_type['right'] == 'reflective':
            im1 = i-1
            ip1 = im1
            u[i] = - u_nm1[i] + 2*u_n[i] + \
                C2*(0.5*(c[i]**2 + c[ip1]**2)*(u_n[ip1] - u_n[i])  - \
                    0.5*(c[i]**2 + c[im1]**2)*(u_n[i] - u_n[im1])) + \
            dt2*f(x[i], t[n])
        elif bc_type['right'] == 'non-reflective':
            im1 = i-1
            ip1 = im1
            u[i] = - u_nm1[i] + 2*u_n[i] + \
                C2*(0.5*(c[i]**2 + c[ip1]**2)*(u_n[ip1] - u_n[i])  - \
                    0.5*(c[i]**2 + c[i]**2)*(u_n[i] - u_nm1[i])) + \
            dt2*f(x[i], t[n])
        elif bc_type['right'] == 'harmonic':
            u[i] = U_L(t[n+1])   
        else :
            raise ValueError('bc right error : {} not implemented'.format(bc_type['right'])) 

        
    return u, u_n, u_nm1