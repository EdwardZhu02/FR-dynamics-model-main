import numpy as np



s_maxlayer =  11                        # Number of soil layers

i_cz = np.float64()                     # Average soil elevation (m)
i_delz = np.float64()                   # Thickness of each soil layer (m)

s_delz = np.empty(shape=(n,), dtype=np.float64)     # Thickness of each soil layer (m)

cH2Ol_s = np.empty(shape=(n,), dtype=np.float64)    # Soil water content in each layer (could also be called theta(:)) (m)

su__ = np.empty(shape=(n,), dtype=np.float64)       # Soil saturation degree in each layer (-)

s_thetar = np.empty(shape=(n,), dtype=np.float64)   # Residual soil water content (-)

s_thetas = np.empty(shape=(n,), dtype=np.float64)   # Saturated soil water content (-)

