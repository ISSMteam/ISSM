import numpy as np


def shock(x, y):

    vel = np.exp(-(np.sqrt((x + 0.1)**2 + (y + 0.1)**2) - 0.75)**2 * 1e6) + ((x + 0.1)**2 + (y + 0.1)**2) / 2
    return vel
