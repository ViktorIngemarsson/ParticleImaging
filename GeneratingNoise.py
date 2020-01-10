import numpy as np
from itertools import product, count

def GenerateUnitVectors(n):
    '''
    Generates matrix NxN of unit length vectors
    '''
    phi = np.random.uniform(0, 2*np.pi, (n, n))
    v = np.stack((np.cos(phi), np.sin(phi)), axis=-1)
    return v

def Qint(t):
    '''
    Quintic interpolation
    '''
    return t * t * t * (t * (t * 6 - 15) + 10)

def Generate2DPerlinNoise(size, ns):
    '''
    generate_2D_perlin_noise(size, ns)
    Generate 2D array of size x size filled with Perlin noise.
    Parameters
    ----------
    size : int
        Size of 2D array size x size.
    ns : int
        Distance between nodes.
    Returns
    -------
    m : ndarray
        The 2D array filled with Perlin noise.
    '''
    nc = int(size / ns)
    grid_size = int(size / ns + 1)
    v = GenerateUnitVectors(grid_size)
    ad, ar = np.arange(ns), np.arange(-ns, 0, 1)
    vd = np.zeros((ns, ns, 4, 1, 2))
    for (l1, l2), c in zip(product((ad, ar), repeat=2), count()):
        vd[:, :, c, 0] = np.stack(np.meshgrid(l2, l1, indexing='xy'), axis=2)

    d = Qint(np.stack((np.zeros((ns, ns, 2)),
                     np.stack(np.meshgrid(ad, ad, indexing='ij'), axis=2)),
           axis=2) / ns)
    d[:, :, 0] = 1 - d[:, :, 1]
    d0 = d[..., 0].copy().reshape(ns, ns, 1, 2)
    d1 = d[..., 1].copy().reshape(ns, ns, 2, 1)
    m = np.zeros((size, size))
    t = m.reshape(nc, ns, nc, ns)

    for i, j in product(np.arange(nc), repeat=2):  # loop through the grid
        av = v[i:i+2, j:j+2].reshape(4, 2, 1)
        at = np.matmul(vd, av).reshape(ns, ns, 2, 2)
        t[i, :, j, :] = np.matmul(np.matmul(d0, at), d1).reshape(ns, ns)
    return m

def FinalNoiseGeneration(resolution, noiseIntensity):
    '''
    Generates background noise using three separate perlin noise generation
    :param resolution: Resolution of quadratic image as pixel on one side
    :param noiseIntensity: Scale of how dark the noise should be at most, 1 is real black.
    :return: Background noise as a nparray
    '''
    img = Generate2DPerlinNoise(resolution, resolution)
    img2 = Generate2DPerlinNoise(resolution, int(resolution/2))
    img3 = Generate2DPerlinNoise(resolution, int(resolution/100))
    rs = np.random.RandomState()
    shot_noise = rs.poisson(10000, (resolution, resolution))
    img = img + img2 + img3 + shot_noise
    img = (img - img.min())/(img.max()-img.min())*noiseIntensity
    return(img)
