import math
from Vector import Vector
import numpy as np

def dotproduct(v1, v2):
    return sum((a * b) for a, b in zip(v1, v2))

def length(v):
    return math.sqrt(dotproduct(v, v))


def angle(v1, v2):
    return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

class SLine:
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2

    def translate(self, ln, dp):
        # TRANSLATE 3D translation of line set
        #
        # LNt = TRANSLATE(LN,dP) translates set of lines LN by dP.
        #   If dP is a Point, the translation corresponds to the
        #   coordinates X, Y and Z.
        #   If dP is a Vector, the translation corresponds to the
        #   components Vx, Vy and Vz.
        #
        # See also SLine, Point, Vector.

        ln_t = ln
        ln_t.p1 = ln_t.p1.translate(dp)
        ln_t.p2 = ln_t.p2.translate(dp)
        return ln_t

    def xrotation(self,ln, phi):
        # XROTATION Rotation around x-axis of line set
        #
        # LNr = XROTATION(LN,phi) rotates set of lines LN around x-axis
        #   by an angle phi [rad].
        #
        # See also SLine.

        ln_r = ln
        ln_r.p1 = ln_r.p1.xrotation(phi)
        ln_r.p2 = ln_r.p2.xrotation(phi)
        return ln_r

    def yrotation(self,ln,phi):
        # YROTATION Rotation around y-axis of line set
        #
        # LNr = YROTATION(LN,phi) rotates set of lines LN around y-axis
        #   by an angle phi [rad].
        #
        # See also SLine.

        ln_r = ln
        ln_r.p1 = ln_r.p1.yrotation(phi)
        ln_r.p2 = ln_r.p2.yrotation(phi)
        return ln_r

    def zrotation(self,ln,phi):
        # ZROTATION Rotation around z-axis of line set
        #
        # LNr = ZROTATION(LN,phi) rotates set of lines LN around z-axis
        #   by an angle phi [rad].
        #
        # See also SLine.

        ln_r = ln
        ln_r.p1 = ln_r.p1.zrotation(phi)
        ln_r.p2 = ln_r.p2.zrotation(phi)
        return ln_r

    def numel(self, ln):
        # NUMEL umber of lines
        #
        # N = NUMEL(LN) number of lines in set LN.
        #
        # See also SLine.

        return ln.p1.size

    def size(self, ln, varargin):
        # SIZE Size of line set
        #
        # S = SIZE(LN) returns a two-element row vector with the number
        #   of rows and columns in the line set LN.
        #
        # S = SIZE(LN,DIM) returns the length of the dimension specified
        #   by the scalar DIM in the line set LN.
        #
        # See also SLine.

        if varargin is None:
            s = ln.p1.size()
        else:
            s = ln.p1.size(varargin[1])
        return s

    def angle(self,ln1,ln2):
        # ANGLE Angle (Components)
        #
        # PHI = ANGLE(LN1,LN2) calculates the angle between the set of
        #   lines LN1 and LN2.
        #
        # See also SLine.

        phi = angle(ln1.p2 - ln1.p1, ln2.p2 - ln2.p1)
        return phi

    def versor(self, ln):
        # VERSOR Unitary vector
        #
        # U = VERSOR(LN) returns the unit vector set corresponding to
        #   the line set LN.
        #   The coordinates X, Y and Z of U are the ones of the initial
        #   point of LN.
        #
        # See also SLine.

        v = Vector(np.zeros(ln.shape()), np.zeros(ln.shape()), np.zeros(ln.shape()), ln.p2.X - ln.p1.X, ln.p2.Y - ln.p1.Y,
               ln.p2.Z - ln.p1.Z)
        u = v.normalize()
        return u

    def tovector(self,ln):
        # TOVECTOR Line to vector
        #
        # V = TOVECTOR(LN) converts the set of lines LN into the set of
        #   vectors V. The coordinates X, Y and Z of V correspond to
        #   the coordinates of the first point of LN and the components
        #   Vx, Vy and Vz of V correspond to the difference between the
        #   coordinates of the first and second point of LN.
        #
        # See also SLine, Vector.

        v = Vector(ln.p1.X, ln.p1.Y, ln.p1.Z, ln.p2.X - ln.p1.X, ln.p2.Y - ln.p1.Y, ln.p2.Z - ln.p1.Z)
        return v