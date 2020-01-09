import numpy as np


def GeneratingCoordinates(rho, sphere_r, x0, y0):

    """
    Generates coordinates in a grid-like fashion within a circle (x,y)
    :param rho: Ray density per square meter.
    :param sphere_r: Particle radius [m].
    :param x0: Coordinate for center of sphere
    :param y0: Coordinate for center of sphere
    :return: Grid-like coordinates within a circle
    """

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

    return x * sphere_r + x0, y * sphere_r + y0
