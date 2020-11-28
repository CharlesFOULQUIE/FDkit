DataSet={
"Probleme":{"equation":"Wave"},                # 'Diffusion', 'DiffAdvec', 'Burger', 'inviscidBurger', 'Advec', 'Wave'

"Domain":{"time": 7.0, 
          "length": 2.0},
          
"Discretization":{"spatial resolution": 600,        # spatial resolution
                  "fourier number":0.5,             # Maximum Fourier number
                  "courant number":0.8,             # Maximum Courant number
                  "scheme":"CD",              # Numerical scheme
                  "version":"scalar"},          # 'scalar' or 'vectorized'
                  
"Animation":{"active":True,                         # Desactivate PlotAndStore
             "skip":5},                             # skip PlotAndStore 
             
"Init":{"active":False, 
            "shape":"gaussian",                     # 'gaussian','cosinehat', 'half-cosinehat', 'plug' or 'sinus'
            "loc":0.5,                              #  location of initial condition
            "sigma":0.05},                          #  width measure of the impulsion
        
"Source":{"active":False,
          "type":"progressive_gaussian_peak",       # 'progressive_gaussian_peak', 'harmonic_gaussian_peak'
          "loc":0.5,                                # location of initial condition
          "sigma":0.05},                            # width measure of the source
          
"Medium":{"diffusion_coeff":0.005,                   # diffusion coeff
          "celerity":1.0,                            # medium velocity
          "active_zone":False,                        # Active zone with slowness factor
          "zone":[1.0, 1.25],                        # interval for right medium
          "slowness_factor":3.0},                    # inverse of wave vel. in right medium

"BoundaryConditions":{"left":'harmonic',           # x =0, 'reflective', 'non-reflective', 'periodic', 'specified'  
                      "right":'non-reflective'},      # x =L, 'reflective', 'non-reflective', 'periodic',  'specified'
}



# Burger FE_LW
# Diffusion / FE_CD
# DiffAdvec / FE_LW
# Advec / Leapfrog, LF, LW, UP, FE_CD 
# Wave / CD 





