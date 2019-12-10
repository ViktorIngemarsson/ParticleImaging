from Spherical import Spherical
import copy
#import numpy as np
import math

class ParticleSpherical:
    def __init__(self, c, r, nm, np):
        self.sp = Spherical(c, r)
        self.nm = nm
        self.np = np

    def translate(self,dp):
        bead_t = copy.deepcopy(self)
        return bead_t.sp.translate(dp)

    def xrotation(self, phi):
        bead_r = copy.deepcopy(self)
        return bead_r.sp.xrotation(phi)

    def yrotation(self, phi):
        bead_r = copy.deepcopy(self)
        return bead_r.sp.yrotation(phi)

    def zrotation(self, phi):
        bead_r = copy.deepcopy(self)
        return bead_r.sp.zrotation(phi)

    def numel(self):
        return self.sp.numel()

    def size(self, varargin = None):
        if varargin is None:
            s = self.sp.size()
        else:
            s = self.sp.size(varargin[1])
        return s

    def barycenter(self):
        # BARYCENTER Spherical particle center of mass
        #
        # P = BARYCENTER(BEAD) returns the point P representing the
        #   center of mass of the spherical particle BEAD.
        #
        # See also ParticleSpherical, Point, Spherical.
        return self.sp.c

    def scattering(self,r,err = None,N = None):
        # SCATTERING Scattered rays
        #
        # S = SCATTERING(BEAD,R) calculates the set of scattered rays S
        #   due to the scattering of the set of rays R on the spherical
        #   particle BEAD.
        #   S is a structure indexed on the scattering events. S(n).r is
        #   the n-th reflected set of rays and S(n).t is the n-th
        #   transmitted set of rays.
        #
        # S = SCATTERING(BEAD,R,ERR) stops the calculation when the
        #   remaining power of the scattered rays is less than ERR
        #   times the power of the incoming rays [default ERR=1e-12].
        #
        # S = SCATTERING(BEAD,R,ERR,N) stops the calculation when the
        #   remaining power of the scattered rays is less than ERR
        #   times the power of the incoming rays [default ERR=1e-12] or
        #   the number of iterations is N [default N=10].
        #
        # See also ParticleSpherical, Ray.
        if (N == None):
            N = 10
        elif (err == None):
            err = 1e-12
        a = r.snellslaw(self.sp, self.nm, self.np, 1)
        c = [a]
        for n in range(0, N-1, 1):
            l = c[n]['r_r'].snellslaw(self.sp, self.np, self.nm, 2)
            c.append(l)
            if c[n + 1]['r_r'].P < r.P * err or math.isnan(c[n + 1]['r_r'].P):
                break
        return c


