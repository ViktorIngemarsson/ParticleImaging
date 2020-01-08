import math
import random
import numpy as np
import matplotlib.pyplot as plt

def GeneratingCoordinates(rho, sphere_r, x0, y0):
    n = int(np.ceil(np.sqrt(rho) * sphere_r))
    m = 0
    xt = np.zeros(n * n)
    yt = np.zeros(n * n)
    xq = np.zeros(n * n)
    yq = np.zeros(n * n)
    for i in range(1, n+1):
        for j in range(1, n+1):
            k = (i-1)*n + j - 1
            xt[k] = 1/n*i
            yt[k] = 1/n*j
            if np.sqrt(np.power(xt[k], 2)+np.power(yt[k], 2)) <= 1:
                xq[m] = xt[k]
                yq[m] = yt[k]
                m += 1

    yvt = np.arange(-1, 0, 1 / n)
    yvb = np.arange(1 / n, 1+1/n, 1 / n)
    xvt = np.zeros(n)
    xvb = np.zeros(n)

    yhl = np.zeros(n)
    yhr = np.zeros(n)
    xhl = np.arange(-1, 0, 1/n)
    xhr = np.arange(1/n, 1+1/n, 1/n)

    xc = np.arange(0, 1, 1)
    yc = np.arange(0, 1, 1)

    x = np.concatenate((xq[0: m], -xq[0: m], -xq[0: m], xq[0: m], xvt, xvb, xhl, xhr, xc))
    y = np.concatenate((yq[0: m], yq[0: m], -yq[0: m], -yq[0: m], yvt, yvb, yhl, yhr, yc))

    plt.scatter(x,y)
    # plt.show()
    return x * sphere_r + x0, y * sphere_r + y0


  #  alpha = 0
  #  b = np.round(alpha * np.sqrt(n))
  #  phi = (np.sqrt(5) + 1) / 2
  #  r = np.empty(n)
  #  theta = np.empty(n)
  #  for k in range(1, n):
  #      r[k] = radius(k, n, b)
  #      theta[k] = 2 * np.pi * k / phi ** 2
  #  maxRadius = max(r)
  #  r = r / maxRadius
  #  x = r * np.cos(theta)
  #  y = r * np.sin(theta)
  #  return x * sphere_r, y * sphere_r


def radius(k, n, b):
    if k > n - b:
        r = 1
    else:
        r = np.sqrt(k - 1 / 2) / np.sqrt(n - (b + 1) / 2)
    return r
