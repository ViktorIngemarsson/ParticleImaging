import numpy as np
from Point import Point
from SLine import SLine
import copy


class Plane:
    def __init__(self, point0, point1, point2):
        self.p0 = point0
        self.p1 = point1
        self.p2 = point2

    def disp(self):
        print(self.p0.X)
        print(self.p0.Y)
        print(self.p0.Z)
        print(self.p1.X)
        print(self.p1.Y)
        print(self.p1.Z)
        print(self.p2.X)
        print(self.p2.Y)
        print(self.p2.Z)

    def translate(self, dp):
        # TRANSLATE 3D translation of plane set
        # PLt = TRANSLATE(PL,dP) translates set of planes PL by dP.
        #   If dP is a Point, the translation corresponds to the
        #   coordinates X, Y and Z.
        #   If dP is a Vector, the translation corresponds to the
        #   components Vx, Vy and Vz.
        ## See also Plane, Point, Vector.
        pl_t = copy.deepcopy(self)
        pl_t.p0 = pl_t.p0.translate(dp)
        pl_t.p1 = pl_t.p1.translate(dp)
        pl_t.p2 = pl_t.p2.translate(dp)
        return pl_t

    def xrotation(self, phi):
        pl_r = copy.deepcopy(self)
        pl_r.p0 = pl_r.p0.xrotation(phi)
        pl_r.p1 = pl_r.p1.xrotation(phi)
        pl_r.p2 = pl_r.p2.xrotation(phi)
        return pl_r

    def yrotation(self, phi):
        pl_r = copy.deepcopy(self)
        pl_r.p0 = pl_r.p0.yrotation(phi)
        pl_r.p1 = pl_r.p1.yrotation(phi)
        pl_r.p2 = pl_r.p2.yrotation(phi)
        return pl_r

    def zrotation(self, phi):
        pl_r = copy.deepcopy(self)
        pl_r.p0 = pl_r.p0.zrotation(phi)
        pl_r.p1 = pl_r.p1.zrotation(phi)
        pl_r.p2 = pl_r.p2.zrotation(phi)
        return pl_r

    def numel(self):
        return self.p0.size()

    def size(self, varargin=None):
        # SIZE Size of the plane set
        #
        # S = SIZE(PL) returns a two-element row vector with the number
        #   of rows and columns in the plane set PL.
        #
        # S = SIZE(PL,DIM) returns the length of the dimension specified
        #   by the scalar DIM in the plane set P .
        #
        # See also Plane.
        if varargin is None:
            s = self.p0.size()
        else:
            s = self.p0.size(varargin[1])
        return s

    def intersectionpoint(self, d):
        # INTERSECTIONPOINT Intersection point between plane and line/vector/ray
        #
        # P = INTERSECTIONPOINT(PL,D) calculates intersection Points
        #   between a set of lines (or vectors) D and the set of planes PL.
        #   If D is parallel to PL, the coordiantes of P are NaN.
        #
        # See also Plane, Point, SLine, Vector, Ray.
        pl = copy.deepcopy(self)

        if isinstance(d, SLine):
            ln = copy.deepcopy(d)
        else:
            ln = d.toline()

        c1 = pl.p1.minus(pl.p0)
        c2 = pl.p2.minus(pl.p0)

        c0 = c1.mtimes(c2)

        c = ln.p2.minus(ln.p1)

        t = c0.times(pl.p0.minus(ln.p1)) / (c0.times(c))  # pl.p0.minus(ln.p1)   (c0.mtimes(c))

        # p = np.add(ln.p1, np.multiply(t, c))
        return ln.p1.plus(c.times(t))

    def perpline(self, p):
        # PERPLINE Line perpendicular to plane passing by point
        #
        # LN = PERPLINE(PL,P) calculates the line set LN perpendicular
        #   to the plane set PL and passing by the point set P.
        #
        # See also Plane.
        pl = copy.deepcopy(self)

        c1 = pl.p1.minus(pl.p0)
        pl = copy.deepcopy(self)
        c2 = pl.p2.minus(pl.p0)

        c0 = c1.mtimes(c2)

        lnp1 = Point(np.multiply(p.X, np.ones(c0.size())), np.multiply(p.Y, np.ones(c0.size())),
                     np.multiply(p.Z, np.ones(c0.size())))

        lnp2 = lnp1.plus(c0)

        ln = SLine(lnp1, lnp2)

        return ln

    def contains(self, p, ln):
        # CONTAINS Plane containing point and line (Static)
        #
        # PL = CONTAINS(P,LN) calculates plane set PL contianing
        #   the point set P and the line set LN.
        #
        # See also Plane, Point, SLine.

        plp0 = Point(np.multiply(p.X, np.ones(ln.shape())), np.multiply(p.Y, np.ones(ln.shape())),
                     np.multiply(p.Z, np.ones(ln.shape())))
        plp1 = Point(np.multiply(ln.p1.X, np.ones(p.shape())), np.multiply(ln.p1.Y, np.ones(p.shape())),
                     np.multiply(ln.p1.Z, np.ones(p.shape())))
        plp2 = Point(np.multiply(ln.p2.X, np.ones(p.shape())), np.multiply(ln.p2.Y, np.ones(p.shape())),
                     np.multiply(ln.p2.Z, np.ones(p.shape())))
        pl = Plane(plp0, plp1, plp2)

        return pl

    def perpto(self, ln, p):
        # PERPTO Plane perpendicular to line and passing by point (Static)
        #
        # PL = PERPTO(LN,P) calculates plane set PL perpendicular
        #   to line set LN and passing by point set P.
        #   LN and P must have the same size, or just must be a singleton.
        #
        # See also Plane, Sline, Point.

        ln.p1.X = np.multiply(ln.p1.X, np.ones(p.shape()))
        ln.p1.Y = np.multiply(ln.p1.Y, np.ones(p.shape()))
        ln.p1.Z = np.multiply(ln.p1.Z, np.ones(p.shape()))
        ln.p2.X = np.multiply(ln.p2.X, np.ones(p.shape()))
        ln.p2.Y = np.multiply(ln.p2.Y, np.ones(p.shape()))
        ln.p2.Z = np.multiply(ln.p2.Z, np.ones(p.shape()))

        p.X = np.multiply(p.X, np.ones(ln.shape()))
        p.Y = np.multiply(p.Y, np.ones(ln.shape()))
        p.Z = np.multiply(p.Z, np.ones(ln.shape()))

        c0 = ln.translate(-ln.p1).p2

        Xt = np.zeros(c0.shape())
        Yt = np.zeros(c0.shape())
        Zt = np.ones(c0.shape())
        Xt[np.power(c0.Z, 2) > np.multiply(100, (np.add(np.power(c0.X, 2), np.power(c0.Y, 2))))] = 1
        Zt[np.power(c0.Z, 2) > np.multiply(100, np.add(np.power(c0.X, 2), np.power(c0.Y, 2)))] = 0
        c1 = c0 * Point(Xt, Yt, Zt)

        c2 = c0 * c1

        pl = Plane(p, p + c1, p + c2)

        return pl
