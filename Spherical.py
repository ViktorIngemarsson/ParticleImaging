import numpy as np
from Point import Point
from SLine import SLine
from Plane import Plane
import copy

class Spherical:
    def __init__(self, centers, radiuses):
        '''
        Generates set of spherical objects
        :param centers: centers of sphere as point objects
        :param radiuses: radius of spheres
        '''
        self.c = centers
        self.r = radiuses

    def disp(self):
        '''
        Displays spherical objects state
        '''
        print(self.c.X)
        print(self.c.Y)
        print(self.c.Z)
        print(self.r)

    def translate(self, dp):
        '''
        Translate position of spherical object
        If dp is point the translation corresponds to the coordinates X, Y and Z.
        If dP is a Vector, the translation corresponds to the components Vx, Vy and Vz.
        :param dp: point or vector to translate set of spherical objects
        :return: new set of spherical objects
        '''
        sp = copy.deepcopy(self)
        sp_t = sp
        sp_t.c = sp.c.translate(dp)

        return sp_t

    def numel(self):
        '''
        Calculates number of spheres in set
        :return: number of spheres [int]
        '''
        sp = copy.deepcopy(self)
        return sp.c.size()

    def intersectionpoint(self, d, n):
        '''
        Calculates intersection points between spheres and line/vector/ray
        If d does not intersect with sphericals, coordinates of return are NaN
        :param d: vectors or lines
        :param n: what intersectionpoint to use, first or second [1,2]
        :return: intersectionpoints
        '''
        sp = copy.deepcopy(self)

        if isinstance(d, SLine):
            ln = d
        else:
            ln = d.toline()

        lnc = ln.p2.minus(ln.p1)

        A = lnc.times(lnc)
        B = ln.p1.minus(sp.c).times(lnc.times(2))
        C = ln.p1.minus(sp.c).times(ln.p1.minus(sp.c)) - np.power(sp.r, 2)

        delta = np.power(B, 2) - np.multiply(4, np.multiply(A, C))

        if n == 1:
            t1 = np.divide(-B - np.sqrt(delta), 2*A)
        else:
            t1 = np.divide(-B + np.sqrt(delta), 2 * A)

        p = ln.p1.plus(lnc.times(t1))

        if (delta < 0).any():
            p.X[delta < 0] = np.nan
            p.Y[delta < 0] = np.nan
            p.Z[delta < 0] = np.nan

        return p

    def perpline(self, p):
        '''
        Calculates the line perpendicular to spheres surface at point p.
        :param p: point p
        :return: line perpendicular to spherical surface
        '''
        sp = copy.deepcopy(self)
        p1 = Point(np.multiply(sp.c.X, np.ones(np.shape(p))), np.multiply(sp.c.Y, np.ones(np.shape(p))), np.multiply(sp.c.Z, np.ones(np.shape(p))))
        p2 = Point(np.multiply(p.X, np.ones(np.shape(p))), np.multiply(p.Y, np.ones(np.shape(p))), np.multiply(p.Z, np.ones(np.shape(p))))
        ln = SLine(p1, p2)
        return ln

    def tangentplane(self, p):
        '''
        Calculates tangentplane between point and spherical object self
        :param p: point
        :return: return plane
        '''
        sp = copy.deepcopy(self)
        pl = Plane(0, 0, 0).perpto(sp.perpline(p), p)
        return pl

