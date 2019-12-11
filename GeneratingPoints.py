import math
import random
import numpy as np

def GeneratingCoordinates(n):
    alpha = 0
    b = np.round(alpha * np.sqrt(n))
    phi = (np.sqrt(5) + 1) / 2
    r = np.empty(n)
    theta = np.empty(n)
    for k in range(1,n):
        r[k] = radius(k, n, b)
        theta[k] = 2 * np.pi * k / phi**2
    maxRadius = max(r)
    r = r/maxRadius
    x = r*np.cos(theta)
    y = r * np.sin(theta)
    return x,y

def radius(k, n, b):
    if k > n - b:
        r = 1
    else:
        r = np.sqrt(k - 1 / 2) / np.sqrt(n - (b + 1) / 2)
    return r