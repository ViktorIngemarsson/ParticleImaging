import Point as Point
from SLine import SLine
import numpy as np

class Vector:
    def __init__(self, xcomponent, ycomponent, zcomponent):
        self.Vx = xcomponent
        self.Vy = ycomponent
        self.Vz = zcomponent

    def vector(self, x, y, z):
        # VECTOR(X, Y, Z, Vx, Vy, Vz) constructs a set of vectorsa
        # with coordinates x, y and z adn components Vx, Vy and Vz.
        # X, Y, Z, Vx, Vy and Vz must be real scalar matrices with the same size.
        #
        # See also Vector, Point.
        v = Point.Point(x, y, z)
        v.Vx = self.Vx
        v.Vy = self.Vy
        v.Vz = self.Vz
        return v

    @staticmethod
    def xrotation(self, v, phi):
        # XROTATION Rotation around x - axis of vector set
        #
        # Vr = XROTATION(V, phi) rotates set of vectors V around x - axis
        # by an angle phi[rad].
        # It rotates both the coordinates X, Y and Z and the components Vx, Vy and Vz.
        # See also Vector.

        v_r = Point(v,phi)
        v_r.Vy = v.Vy*np.cos(phi) - v.Vz*np.sin(phi)
        v_r.Vz = v.Vy*np.sin(phi) + v.Vz*np.cos(phi)
        return v_r

    @staticmethod
    def yrotation(self, v, phi):
        # YROTATION Rotation around x - axis of vector set
        #
        # Vr = YROTATION(V, phi) rotates set of vectors V around y - axis
        # by an angle phi[rad].
        # It rotates both the coordinates X, Y and Z and the components Vx, Vy and Vz.
        # See also Vector.

        v_r = Point(v, phi)
        v_r.Vx = v.Vx * np.cos(phi) + v.Vz * np.sin(phi)
        v_r.Vz = -v.Vx * np.sin(phi) + v.Vz * np.cos(phi)
        return v_r

    @staticmethod
    def zrotation(self, v, phi):
        # ZROTATION Rotation around x - axis of vector set
        #
        # Vr = ZROTATION(V, phi) rotates set of vectors V around z - axis
        # by an angle phi[rad].
        # It rotates both the coordinates X, Y and Z and the components Vx, Vy and Vz.
        # See also Vector.

        v_r = Point(v, phi)
        v_r.Vx = v.Vx * np.cos(phi) - v.Vy * np.sin(phi)
        v_r.Vy = v.Vx * np.sin(phi) + v.Vy * np.cos(phi)

        return v_r

    @staticmethod
    def uminus(v):
        # UMINUS Unitary minus(components)
        #
        # Vm = UMINUS(V) Unitary minus.
        # -V inverts the components Vx, Vy and Vz of V.
        # The coordinates X, Y and Z are left unchanged.
        #
        # See also Vector.

        v_m = v
        v_m.Vx = -v_m.Vx
        v_m.Vy = -v_m.Vy
        v_m.Vz = -v_m.Vz

        return v_m

    @staticmethod
    def plus(self, v1, v2):
        # PLUS Binary addition(components)
        #
        # V = PLUS(V1, V2) Binary addition(V=V1 + V2).
        # The components Vx, Vy and Vz of V are the sum of the
        # components of V1 and V2.
        # The coordinates X, Y and Z of V are the ones of V1.
        #
        # See also Vector.

        v = v1
        v.Vx = v1.Vx + v2.Vx
        v.Vy = v1.Vy + v2.Vy
        v.Vz = v1.Vz + v2.Vz

        return v

    @staticmethod
    def minus(self, v1, v2):
        # PLUS Binary addition(components)
        #
        # V = MINUS(V1, V2) Binary subtraction(V=V1 - V2).
        # The components Vx, Vy and Vz of V are the difference of the
        # components of V1 and V2.
        # The coordinates X, Y and Z of V are the ones of V1.
        #
        # See also Vector.

        v = v1
        v.Vx = v1.Vx - v2.Vx
        v.Vy = v1.Vy - v2.Vy
        v.Vz = v1.Vz - v2.Vz

        return v

    @staticmethod
    def mtimes(self, a, b):
        # MTIMES Vector product(components)
        #
        # V = MTIMES(V1, V2) Vector product(V=V1 * V2).
        # V is a Vector whose components Vx, Vy and Vz are the vector
        # product of the components of V1 and V2.
        # The coordinates X, Y and Z of V are the ones of V1.
        #
        # V = MTIMES(A, V2) Product by scalar(P=A * V2).
        # V is a Vector whose components Vx, Vy and Vz are the components
        # of V2 multiplied by the scalar( or scalar matrix) A.
        # The coordinates X, Y and Z of V are the ones of V2.
        #
        # V = MTIMES(P1, B) Product by scalar(P=P1 * B).
        # V is a Vector whose components Vx, Vy and Vz are the components
        # of V1 multiplied by the scalar( or scalar matrix) B.
        # The coordinates X, Y and Z of V are the ones of V1.
        #
        # See also Vector.

        if a.isinstance(Vector) & b.isinstance(Vector):
            v1 = a
            v2 = b
            m = v1
            m.Vx = np.dot(v1.Vy, v2.Vz)-np.dot(v1.Vz, v2.Vy)
            m.Vy = -np.dot(v1.Vx, v2.Vz)+np.dot(v1.Vz, v2.Vx)
            m.Vz = np.dot(v1.Vx, v2.Vy)-np.dot(v1.Vy, v2.Vx)
        elif a.isinstance(Vector):
            v1 = a
            m = v1
            m.X = m.X * np.ones([b.shape], dtype = int)
            m.Y = m.Y * np.ones([b.shape], dtype = int)
            m.Z = m.Z * np.ones([b.shape], dtype = int)
            m.Vx = np.dot(v1.Vx, b)
            m.Vy = np.dot(v1.Vy, b)
            m.Vz = np.dot(v1.Vz, b)
        elif b.isinstance(Vector):
            v2 = b
            m = v2
            m.X = m.X * np.ones([a.shape], dtype = int)
            m.Y = m.Y * np.ones([a.shape], dtype = int)
            m.Z = m.Z * np.ones([a.shape], dtype = int)
            m.Vx = np.dot(a, v2.Vx)
            m.Vy = np.dot(a, v2.Vy)
            m.Vz = np.dot(a, v2.Vz)
        else:
            m = np.dot(a, b)

        return m

    @staticmethod
    def times(self, a, b):
        # TIMES Scalar product(components)
        #
        # M = TIMES(V1, V2) Scalar product(M=V1. * V2).
        # M is a scalar matrix obtained by the scalar product
        # of the components Vx, Vy and Vz of V1 and V2.
        #
        # V = TIMES(A, V2) Product by scalar(V=A. * V2).
        # V is a Vector whose components Vx, Vy and Vz are the components
        # of V2 multiplied by the scalar( or scalar matrix) A.
        # The coordinates X, Y and Z of V are the ones of V2.
        #
        # V = TIMES(V1, B) Product by scalar(V=V1. * B).
        # V is a Vector whose components Vx, Vy and Vz are the components
        # of V1 multiplied by the scalar( or scalar matrix) B.
        # The coordinates X, Y and Z of V are the ones of V1.
        #
        # See also Vector.

        if a.isinstance(Vector) & b.isinstance(Vector):
            v1 = a
            v2 = b
            m = np.dot(v1.Vx, v2.Vx) + np.dot(v1.Vy, v2.Vy) + np.dot(v1.Vz, v2.Vz)
        elif a.isinstance(Vector):
            v1 = a
            m = v1
            m.X = m.X * np.ones([b.shape], dtype=int)
            m.Y = m.Y * np.ones([b.shape], dtype=int)
            m.Z = m.Z * np.ones([b.shape], dtype=int)
            m.Vx = np.dot(v1.Vx, b)
            m.Vy = np.dot(v1.Vy, b)
            m.Vz = np.dot(v1.Vz, b)
        elif b.isinstance(Vector):
            v2 = b
            m = v2
            m.X = m.X * np.ones([a.shape], dtype=int)  # om fel kanske detta ska göras med np.dot
            m.Y = m.Y * np.ones([a.shape], dtype=int)
            m.Z = m.Z * np.ones([a.shape], dtype=int)
            m.Vx = np.dot(a, v2.Vx)
            m.Vy = np.dot(a, v2.Vy)
            m.Vz = np.dot(a, v2.Vz)
        else:
            m = np.dot(a, b)

        return m

    @staticmethod
    def rdivide(self, v, b):
        # RDIVIDE Right division(components)
        #
        # Vd = RDIVIDE(V, B) Right division(Vd=V. / B).
        # V is a Vector whose components Vx, Vy and Vz are the components
        # of V divided by the scalar( or scalar matrix) B.
        #
        # See also Vector.

        v_d = v
        v_d.Vx = np.divide(v.Vx, b)
        v_d.Vy = np.divide(v.Vy, b)
        v_d.Vz = np.divide(v.Vz, b)

        return v_d

    @staticmethod
    def versor(self, v):
        # VERSOR Unitary vector
        #
        # U = VERSOR(V) returns the unit vector set corresponding to
        # the vector set V.
        # The coordinates X, Y and Z of U are the ones of V.
        #
        # See also Vector.

        u = v.normalize()  # tror inte detta kommer funkar för jag vet inte vad den gör...

        return u

    @staticmethod
    def topoint(self, v):
        # TOPOINT Vector to point
        #
        # P = TOPOINT(V) converts the set of vectors V into the set of
        # points P.The coordinates X, Y and Z of the points are the
        # components of V.
        #
        # See also Vector, Point.

        p = Point(np.real(v.Vx), np.real(v.Vy), np.real(v.Vz))

        return p

    @staticmethod
    def toline(self, v):
        # TOLINE Vector to line
        #
        # LN = TOLINE(V) converts the set of vectors V into the set of
        # lines LN.The coordinates X, Y and Z of the initial points
        # of the lines are the coordinates of V and the coordinates
        # of the final points are the sum of the coordinates of V and
        # of the components of V.
        #
        # See also Vector, SLine.

        ln = SLine(Point(v.X, v.Y, v.Z), Point(v.X+np.real(v.Vx), v.Y+np.real(v.Vy), v.Z+np.real(v.Vz)))

        return ln
