import numpy as np
from Point import Point
from SLine import SLine
import copy


class Plane:
    def __init__(self, point0, point1, point2):
        '''
        Generates a set of planes from three sets of points
        :param point0: a point
        :param point1: a point
        :param point2: a point
        '''
        self.p0 = point0
        self.p1 = point1
        self.p2 = point2

    def disp(self):
        '''
        Displays current state of the set of planes
        '''
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
        '''
        translates set of planes PL by dP.
        If dP is a Point, the translation corresponds to the...
        coordinates X, Y and Z.
        If dP is a Vector, the translation corresponds to the...
        components Vx, Vy and Vz.
        :param dp: point or vector
        :return: translated set of planes
        '''
        pl_t = copy.deepcopy(self)
        pl_t.p0 = pl_t.p0.translate(dp)
        pl_t.p1 = pl_t.p1.translate(dp)
        pl_t.p2 = pl_t.p2.translate(dp)
        return pl_t

    def xrotation(self, phi):
        '''
        Rotates set of planes around x - axis
        :param phi: angle phi[rad] rotated around x - axis
        :return: new set of planes
        '''
        pl_r = copy.deepcopy(self)
        pl_r.p0 = pl_r.p0.xrotation(phi)
        pl_r.p1 = pl_r.p1.xrotation(phi)
        pl_r.p2 = pl_r.p2.xrotation(phi)
        return pl_r

    def yrotation(self, phi):
        '''
        Rotates set of planes around y - axis
        :param phi: angle phi[rad] rotated around y - axis
        :return: new set of planes
        '''
        pl_r = copy.deepcopy(self)
        pl_r.p0 = pl_r.p0.yrotation(phi)
        pl_r.p1 = pl_r.p1.yrotation(phi)
        pl_r.p2 = pl_r.p2.yrotation(phi)
        return pl_r

    def zrotation(self, phi):
        '''
        Rotates set of planes around z - axis
        :param phi: angle phi[rad] rotated around z - axis
        :return: new set of planes
        '''
        pl_r = copy.deepcopy(self)
        pl_r.p0 = pl_r.p0.zrotation(phi)
        pl_r.p1 = pl_r.p1.zrotation(phi)
        pl_r.p2 = pl_r.p2.zrotation(phi)
        return pl_r

    def numel(self):
        '''
        Returns the number of planes in set self.
        :return: Number of planes in set
        '''
        return self.p0.size()

    def size(self, varargin=None):
        '''
        Returns a two-element row vector with the number...
        of rows and columns in the plane set self.
        If given argument returns the length of the dimension specified...
        by the scalar DIM in the plane set self.
        :param varargin: integer of which dim
        :return:length of specified dimension
        '''
        if varargin is None:
            s = self.p0.size()
        else:
            s = self.p0.size(varargin[1])
        return s

    def intersectionpoint(self, d):
        '''
        Calculates intersection Points between a set of lines...
        (or vectors) d and the set of planes self.
        If d is parallel to self, the coordiantes of self are NaN.
        :param d: lines or vectors
        :return: intersecting points
        '''
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
        return ln.p1.plus(c.times(t))

    def perpline(self, p):
        '''
        Calculates the line set ln perpendicular to the...
        plane set self and passing by the point set p.
        :param p: set of points
        :return ln: line set
        '''
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
        '''
        Calculates plane set PL containing..
        the point set p and the line set ln.
        :param p: point set
        :param ln: line set
        :return: new plane set
        '''
        plp0 = Point(np.multiply(p.X, np.ones(ln.size())), np.multiply(p.Y, np.ones(ln.size())),
                     np.multiply(p.Z, np.ones(ln.size())))
        plp1 = Point(np.multiply(ln.p1.X, np.ones(p.size())), np.multiply(ln.p1.Y, np.ones(p.size())),
                     np.multiply(ln.p1.Z, np.ones(p.size())))
        plp2 = Point(np.multiply(ln.p2.X, np.ones(p.size())), np.multiply(ln.p2.Y, np.ones(p.size())),
                     np.multiply(ln.p2.Z, np.ones(p.size())))
        pl = Plane(plp0, plp1, plp2)
        return pl

    def perpto(self, p):
        '''
        Calculates plane set pl perpendicular
        to line set LN and passing by point set P.
        LN and P must have the same size, or just must be a singleton.
        :param p: point set
        :return pl: new plane set
        '''

        ln = copy.deepcopy(self)

        ln.p1.X = np.multiply(ln.p1.X, np.ones(p.size()))
        ln.p1.Y = np.multiply(ln.p1.Y, np.ones(p.size()))
        ln.p1.Z = np.multiply(ln.p1.Z, np.ones(p.size()))
        ln.p2.X = np.multiply(ln.p2.X, np.ones(p.size()))
        ln.p2.Y = np.multiply(ln.p2.Y, np.ones(p.size()))
        ln.p2.Z = np.multiply(ln.p2.Z, np.ones(p.size()))

        p.X = np.multiply(p.X, np.ones(ln.size()))
        p.Y = np.multiply(p.Y, np.ones(ln.size()))
        p.Z = np.multiply(p.Z, np.ones(ln.size()))

        c0 = ln.translate(ln.p1.uminus).p2

        Xt = np.zeros(c0.size())
        Yt = np.zeros(c0.size())
        Zt = np.ones(c0.size())
        #Xt[c0.Z.times(c0.Z) > np.multiply(100, (c0.X.times(c0.X).plus(c0.Y.times(c0.Y))))] = 1
        Xt[np.power(c0.Z, 2) > np.multiply(100, (np.add(np.power(c0.X, 2), np.power(c0.Y, 2))))] = 1
        Zt[np.power(c0.Z, 2) > np.multiply(100, (np.add(np.power(c0.X, 2), np.power(c0.Y, 2))))] = 0
        c1 = c0.mtimes(Point(Xt, Yt, Zt))

        c2 = c0.mtimes(c1)

        pl = Plane(p, p.plus(c1), p.plus(c2))

        return pl
