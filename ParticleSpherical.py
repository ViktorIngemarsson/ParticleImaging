from Spherical import Spherical
import copy
import numpy as np

class ParticleSpherical:
    def __init__(self, c, r, nm, np):
        """
        The constructor constructs a spherical particle with center at point C, radius R, medium refractive index
            nm and particle refractive index np.
        :param c: sphere center coordinates
        :param r: sphere radius [m]
        :param nm: medium refractive index
        :param np: particle refractive index
        """
        self.sp = Spherical(c, r)
        self.nm = nm
        self.np = np

    def scattering(self,r,err = None,N = None):
        """
        scattering() calculates the scattering of the set of rays that interact with the ParticleSpherical object.
        :param r: sphere radius [m]
        :param err: Stops the calculation when the remaining power of the scattered rays is less than err.
        :param N: Stops the calculation when the rays have scattered N times.
        :return: Returns the set of scattered rays.
        """

        if (N == None):
            N = 10
        if (err == None):
            err = 1e-12

        r_r, r_t, perp = r.snellslaw(self.sp, self.nm, self.np, 1)
        rVec = []
        rVec.append({
            "r": r_r,
            "t": r_t
        })
        r_r, r_t, perp = rVec[0]["t"].snellslaw(self.sp, self.np, self.nm, 2)
        rVec.append({
            "r": r_r,
            "t": r_t
        })
        for n in range(1, N, 1):
            r_r, r_t, perp = rVec[n]['r'].snellslaw(self.sp, self.np, self.nm, 2)
            rVec.append({
                "r": r_r,
                "t": r_t
            })
            if (rVec[n + 1]['r'].P < r.P * err).any() or np.isnan(rVec[n + 1]['r'].P).any():
                break

        return rVec


