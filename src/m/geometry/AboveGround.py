import numpy as np

def AboveGround(lat, long, r, height):  # {{{
    r = r + height
    x = r * np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(long))
    y = r * np.cos(np.deg2rad(lat)) * np.sin(np.deg2rad(long))
    z = r * np.sin(np.deg2rad(lat))
# }}}
