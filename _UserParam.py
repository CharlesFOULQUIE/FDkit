DataSet={
"Probleme":{"equation":"Advec"},                # 'Diffusion', 'DiffAdvec', 'Burger', 'inviscidBurger', 'Advec', 'Wave'

"Domain":{"time": 1.0, 
          "length": 1.0},
          
"Discretization":{"spatial resolution": 100,        # spatial resolution
                  "fourier number":0.5,             # Maximum Fourier number
                  "courant number":0.8,             # Maximum Courant number
                  "scheme":"LF",              # Numerical scheme
                  "version":"vectorized"},          # 'scalar' or 'vectorized'
                  
"Animation":{"active":True,                         # Desactivate PlotAndStore
             "skip":5},                             # skip PlotAndStore 
             
"Init":{"active":True, 
            "shape":"gaussian",                     # 'gaussian','cosinehat', 'half-cosinehat', 'plug' or 'sinus'
            "loc":0.5,                              #  location of initial condition
            "sigma":0.05},                          #  width measure of the impulsion
        
"Source":{"active":False,
          "type":"progressive_gaussian_peak",       # 'progressive_gaussian_peak', 'harmonic_gaussian_peak'
          "loc":0.5,                                # location of initial condition
          "sigma":0.05},                            # width measure of the source
          
"Medium":{"diffusion_coeff":0.005,                   # diffusion coeff
          "celerity":1.0,                            # medium velocity
          "active_zone":False,                       # Active zone with slowness factor
          "zone":[1.0, 1.25],                        # interval for right medium
          "slowness_factor":3.0},                    # inverse of wave vel. in right medium

"BoundaryConditions":{"left":'periodic',           # x =0, 'reflective', 'non-reflective', 'periodic', 'specified'  
                      "right":'periodic'},      # x =L, 'reflective', 'non-reflective', 'periodic',  'specified'
}



# Burger FE_LW
# Diffusion / FE_CD
# DiffAdvec / FE_LW
# Advec / Leapfrog, LF, LW, UP, FE_CD 
# Wave / CD 





