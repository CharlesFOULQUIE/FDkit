DataSet={
"Probleme":{"equation":"Diffusion"},                # 'Diffusion', 'DiffAdvec', 'Burger', 'inviscidBurger', 'Advec', 'Wave'

"Domain":{"time": 1.0, 
          "length": 1.0},
          
"Discretization":{"spatial resolution": 100,        # spatial resolution
                  "fourier number":0.5,             # Maximum Fourier number
                  "courant number":1.0,             # Maximum Courant number
                  "scheme":"FE_CD",                 # Numerical scheme
                  "version":"vectorized"},          # 'scalar' or 'vectorized'
                  
"Animation":{"active":True,                         # Desactivate PlotAndStore
             "skip":1},                            # skip PlotAndStore 
             
"Init":{"active":True, 
            "shape":"gaussian",                     # 'gaussian','cosinehat', 'half-cosinehat', 'plug' or 'sinus'
            "loc":0.5,                              #  location of initial condition
            "sigma":0.05},                          #  width measure of the impulsion
        
"Source":{"active":False,
          "type":"progressive_gaussian_peak",       # 'progressive_gaussian_peak', 'harmonic_gaussian_peak'
          "loc":0.5,                                # location of initial condition
          "sigma":0.05},                            # width measure of the source
          
"Medium":{"diffusion_coeff":0.005,                   # diffusion coeff
          "celerity":1.0,                           # medium velocity
          "active_zone":False,                      # Active zone with slowness factor
          "zone":[0.7, 0.8],                        # interval for right medium
          "slowness_factor":4.0}                    # inverse of wave vel. in right medium
}

# Burger FE_LW
# Diffusion / FE_CD
# DiffAdvec / FE_LW
# Advec / Leapfrog, LF, LW, UP, FE_CD 
# Wave / CD 





