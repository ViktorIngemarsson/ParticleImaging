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

    def xrotation(self, phi):
        # XROTATION Rotation around x-axis of sphere set
        #
        # SPr = XROTATION(SP,phi) rotates set of spheres SP around x-axis
        #      by an angle phi [rad]
        #
        # See also Spherical.

        sp = copy.deepcopy(self)
        sp_r = sp
        sp_r.c = sp.c.xrotation(phi)

        return sp_r

    def yrotation(self, phi):
        # YROTATION Rotation around x-axis of sphere set
        #
        # SPr = YROTATION(SP,phi) rotates set of spheres SP around y-axis
        #      by an angle phi [rad]
        #
        # See also Spherical.

        sp = copy.deepcopy(self)
        sp_r = sp
        sp_r.c = sp.c.yrotation(phi)

        return sp_r

    def zrotation(self, phi):
        # ZROTATION Rotation around x-axis of sphere set
        #
        # SPr = ZROTATION(SP,phi) rotates set of spheres SP around z-axis
        #      by an angle phi [rad]
        #
        # See also Spherical.

        sp = copy.deepcopy(self)
        sp_r = sp
        sp_r.c = sp.c.zrotation(phi)

        return sp_r

    def numel(self):
        # NUMEL Number of spheres
        #
        # N = NUMEL(SP) number of spheres in set SP.
        #
        # See also Spherical.

        sp = copy.deepcopy(self)
        return sp.c.size()

    def size(self, varargin):
        # SIZE Size of the sphere set
        #
        # S = SIZE(SP) returns a two-element row vector with the number
        #   of rows and columns in the sphere set PL.
        #
        # S = SIZE(SP,DIM) returns the length of the dimension specified
        #   by the scalar DIM in the sphere set P .
        #
        # See also Spherical.

        sp = copy.deepcopy(self)

        if varargin != None:
            s = np.shape(np.asarray(self.c), varargin[0])
            return s
        else:
            s = np.shape(np.asarray(self.c))
            return s


    def intersectionpoint(self, d, n):
        # INTERSECTIONPOINT Intersection point between sphere and line/vector/ray
        #
        #  P = INTERSECTIONPOINT(SP,D,N) calculates intersection points
        #   between a set of lines (or vectors) D and the set of spheres SP.
        #   The intersection point is selected by  N = {1,2}.
        #   If D does not intersect SP, the coordinates of P are NaN.
        #
        # See also  Spherical, Point, Vector, SLine, Ray.

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
        # PERPLINE Line perpendicular to sphere passing by point
        #
        # LN = PERPLINE(SP,P) calculates the line set LN perpendicular
        #   to the sphere set SP and passing by the point set P.
        #
        # See also Spherical, Point, SLine

        sp = copy.deepcopy(self)

        p1 = Point(np.multiply(sp.c.X, np.ones(np.shape(p))), np.multiply(sp.c.Y, np.ones(np.shape(p))), np.multiply(sp.c.Z, np.ones(np.shape(p))))

        p2 = Point(np.multiply(p.X, np.ones(np.shape(p))), np.multiply(p.Y, np.ones(np.shape(p))), np.multiply(p.Z, np.ones(np.shape(p))))

        ln = SLine(p1, p2)

        return ln

    def tangentplane(self, p):

        sp = copy.deepcopy(self)
        pl = Plane(0, 0, 0).perpto(sp.perpline(p), p)

        return pl

