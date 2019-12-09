import numpy as np
from Point import Point
from SLine import SLine


class Spherical:
    def __init__(self, centers, radiuses):
        self.c = centers
        self.r = radiuses

    def spherical(self, obj=None):
        # SPHERICAL(c, r) constructs a set of spheres
        # with centers c and radii r.
        # c must be a Point and r a real scalar matrix with the same size.
        #
        # See also Spherical, Point.

        obj.c = self.c
        obj.r = self.r

        return obj

    def translate(self, sp, dp):
        # TRANSLATE 3D translation of sphere set
        # SPt = TRANSLATE(SP,dP) translates set of spheres SP by dP
        #    If dP is a Point, the translation corresponds to the
        #   coordinates X, Y and Z.
        #   If dP is a Vector, the translation corresponds to the
        #   components Vx, Vy and Vz.
        #
        # See also Plane, Point, Vector.

        sp_t = sp
        sp_t.c = sp.c.translate(dp)
        return sp_t

    def xrotation(self, sp, phi):
        # XROTATION Rotation around x-axis of sphere set
        #
        # SPr = XROTATION(SP,phi) rotates set of spheres SP around x-axis
        #      by an angle phi [rad]
        #
        # See also Spherical.

        sp_r = sp
        sp_r.c = sp.c.xrotation(phi)

    def yrotation(self, sp, phi):
        # YROTATION Rotation around x-axis of sphere set
        #
        # SPr = YROTATION(SP,phi) rotates set of spheres SP around y-axis
        #      by an angle phi [rad]
        #
        # See also Spherical.

        sp_r = sp
        sp_r.c = sp.c.yrotation(phi)

    def zrotation(self, sp, phi):
        # ZROTATION Rotation around x-axis of sphere set
        #
        # SPr = ZROTATION(SP,phi) rotates set of spheres SP around z-axis
        #      by an angle phi [rad]
        #
        # See also Spherical.

        sp_r = sp
        sp_r.c = sp.c.zrotation(phi)

    def numel(self, sp):
        # NUMEL Number of spheres
        #
        # N = NUMEL(SP) number of spheres in set SP.
        #
        # See also Spherical.

        return sp.c.size

    def size(self, sp, varargin):
        # SIZE Size of the sphere set
        #
        # S = SIZE(SP) returns a two-element row vector with the number
        #   of rows and columns in the sphere set PL.
        #
        # S = SIZE(SP,DIM) returns the length of the dimension specified
        #   by the scalar DIM in the sphere set P .
        #
        # See also Spherical.

        if varargin is None:
            s = sp.c.size()
        else:
            s = sp.c.size(varargin[1])
        return s

    def intersectionpoint(self, sp, d, n):
        # INTERSECTIONPOINT Intersection point between sphere and line/vector/ray
        #
        #  P = INTERSECTIONPOINT(SP,D,N) calculates intersection points
        #   between a set of lines (or vectors) D and the set of spheres SP.
        #   The intersection point is selected by  N = {1,2}.
        #   If D does not intersect SP, the coordinates of P are NaN.
        #
        # See also  Spherical, Point, Vector, SLine, Ray.

        if d.isinstance(SLine):
            ln = d
        else:
            ln = d.toline()

        lnc = ln.p2 - ln.p1

        A = np.multiply(lnc, lnc)
        B = np.multiply(np.multiply(2, ln.p1 - sp.c), lnc)
        C = np.multiply(ln.p1 - sp.c, ln.p1 - sp.c) - np.power(sp.r, 2)

        delta = np.power(B, 2) - np.multiply(4, np.multiply(A, C))

        if n == 1:
            t1 = np.divide(-B - np.sqrt(delta), 2 * A)
        else:
            t1 = np.divide(-B + np.sqrt(delta), 2 * A)

        p = ln.p1 + np.multiply(t1, lnc)

        if delta < 0:
            p.X = np.nan
            p.Y = np.nan
            p.Z = np.nan

        return p

    def perpline(self, sp, p):
        # PERPLINE Line perpendicular to sphere passing by point
        #
        # LN = PERPLINE(SP,P) calculates the line set LN perpendicular
        #   to the sphere set SP and passing by the point set P.
        #
        # See also Spherical, Point, SLine

        p1 = Point(np.multiply(sp.c.X, np.ones(p)), np.multiply(sp.c.Y, np.ones(p)), np.multiply(sp.c.Z, np.ones(p)))

        p2 = Point(np.multiply(p.X, np.ones(p)), np.multiply(p.Y, np.ones(p)), np.multiply(p.Z, np.ones(p)))

        ln = SLine(p1, p2)

        return ln

    def tangentplane(self, sp, p):

        pl = Plane.perpto(sp.perpline(p), p)

        return pl
