from Point import Point
import numpy as np
import copy


class Vector(Point):  # Inherits from Superclass Point
    def __init__(self, X, Y, Z, xcomponent, ycomponent, zcomponent):
        super().__init__(X, Y, Z)
        self.X = X
        self.Y = Y
        self.Z = Z
        self.Vx = xcomponent
        self.Vy = ycomponent
        self.Vz = zcomponent

    def disp(self):
        print(self.X)
        print(self.Y)
        print(self.Z)
        print(self.Vx)
        print(self.Vy)
        print(self.Vz)

    def xrotation(self, phi):
        # XROTATION Rotation around x - axis of vector set
        #
        # Vr = XROTATION(V, phi) rotates set of vectors V around x - axis
        # by an angle phi[rad].
        # It rotates both the coordinates X, Y and Z and the components Vx, Vy and Vz.
        # See also Vector.
        v = copy.deepcopy(self)

        v_r = Point.xrotation(v, phi)
        v_r.Vy = v.Vy * np.cos(phi) - v.Vz * np.sin(phi)
        v_r.Vz = v.Vy * np.sin(phi) + v.Vz * np.cos(phi)
        return v_r

    def yrotation(self, phi):
        # YROTATION Rotation around x - axis of vector set
        #
        # Vr = YROTATION(V, phi) rotates set of vectors V around y - axis
        # by an angle phi[rad].
        # It rotates both the coordinates X, Y and Z and the components Vx, Vy and Vz.
        # See also Vector.
        v = copy.deepcopy(self)

        v_r = Point.yrotation(v, phi)
        v_r.Vx = v.Vx * np.cos(phi) + v.Vz * np.sin(phi)
        v_r.Vz = -v.Vx * np.sin(phi) + v.Vz * np.cos(phi)
        return v_r

    def zrotation(self, phi):
        # ZROTATION Rotation around x - axis of vector set
        #
        # Vr = ZROTATION(V, phi) rotates set of vectors V around z - axis
        # by an angle phi[rad].
        # It rotates both the coordinates X, Y and Z and the components Vx, Vy and Vz.
        # See also Vector.
        v = copy.deepcopy(self)

        v_r = Point.zrotation(v, phi)
        v_r.Vx = v.Vx * np.cos(phi) - v.Vy * np.sin(phi)
        v_r.Vy = v.Vx * np.sin(phi) + v.Vy * np.cos(phi)

        return v_r

    def uminus(self):
        # UMINUS Unitary minus(components)
        #
        # Vm = UMINUS(V) Unitary minus.
        # -V inverts the components Vx, Vy and Vz of V.
        # The coordinates X, Y and Z are left unchanged.
        #
        # See also Vector.

        v_m = copy.deepcopy(self)
        v_m.Vx = -v_m.Vx
        v_m.Vy = -v_m.Vy
        v_m.Vz = -v_m.Vz

        return v_m

    def plus(self, v2):
        # PLUS Binary addition(components)
        #
        # V = PLUS(V1, V2) Binary addition(V=V1 + V2).
        # The components Vx, Vy and Vz of V are the sum of the
        # components of V1 and V2.
        # The coordinates X, Y and Z of V are the ones of V1.
        #
        # See also Vector.

        v = copy.deepcopy(self)
        v.Vx = self.Vx + v2.Vx
        v.Vy = self.Vy + v2.Vy
        v.Vz = self.Vz + v2.Vz

        return v

    def minus(self, v2):
        # PLUS Binary addition(components)
        #
        # V = MINUS(V1, V2) Binary subtraction(V=V1 - V2).
        # The components Vx, Vy and Vz of V are the difference of the
        # components of V1 and V2.
        # The coordinates X, Y and Z of V are the ones of V1.
        #
        # See also Vector.

        v = copy.deepcopy(self)
        v.Vx = self.Vx - v2.Vx
        v.Vy = self.Vy - v2.Vy
        v.Vz = self.Vz - v2.Vz

        return v

    def mtimes(self, b):
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
        a = copy.deepcopy(self)

        if isinstance(a, Vector) & isinstance(b, Vector):
            v1 = copy.deepcopy(a)
            v2 = b
            m = copy.deepcopy(v1)
            m.Vx = np.dot(v1.Vy, v2.Vz) - np.dot(v1.Vz, v2.Vy)
            m.Vy = -np.dot(v1.Vx, v2.Vz) + np.dot(v1.Vz, v2.Vx)
            m.Vz = np.dot(v1.Vx, v2.Vy) - np.dot(v1.Vy, v2.Vx)
        elif isinstance(a, Vector):
            v1 = copy.deepcopy(a)
            m = copy.deepcopy(v1)
            m.X = m.X * np.ones([b.size()], dtype=int)
            m.Y = m.Y * np.ones([b.size()], dtype=int)
            m.Z = m.Z * np.ones([b.size()], dtype=int)
            m.Vx = np.dot(v1.Vx, b)
            m.Vy = np.dot(v1.Vy, b)
            m.Vz = np.dot(v1.Vz, b)
        elif isinstance(b, Vector):
            v2 = copy.deepcopy(b)
            m = copy.deepcopy(v2)
            m.X = m.X * np.ones([a.size()], dtype=int)
            m.Y = m.Y * np.ones([a.size()], dtype=int)
            m.Z = m.Z * np.ones([a.size()], dtype=int)
            m.Vx = np.dot(a, v2.Vx)
            m.Vy = np.dot(a, v2.Vy)
            m.Vz = np.dot(a, v2.Vz)
        else:
            m = np.dot(a, b)

        return m

    def times(self, b):
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
        a = copy.deepcopy(self)

        if isinstance(a, Vector) & isinstance(b, Vector):
            v1 = copy.deepcopy(a)
            v2 = copy.deepcopy(b)
            m = np.dot(v1.Vx, v2.Vx) + np.dot(v1.Vy, v2.Vy) + np.dot(v1.Vz, v2.Vz)
        elif isinstance(a, Vector):
            v1 = copy.deepcopy(a)
            m = copy.deepcopy(v1)
            m.X = m.X * np.ones([b.size()], dtype=int)
            m.Y = m.Y * np.ones([b.size()], dtype=int)
            m.Z = m.Z * np.ones([b.size()], dtype=int)
            m.Vx = np.dot(v1.Vx, b)
            m.Vy = np.dot(v1.Vy, b)
            m.Vz = np.dot(v1.Vz, b)
        elif isinstance(b, Vector):
            v2 = copy.deepcopy(b)
            m = copy.deepcopy(v2)
            m.X = m.X * np.ones([a.size()], dtype=int)  # om fel kanske detta ska göras med np.dot
            m.Y = m.Y * np.ones([a.size()], dtype=int)
            m.Z = m.Z * np.ones([a.size()], dtype=int)
            m.Vx = np.dot(a, v2.Vx)
            m.Vy = np.dot(a, v2.Vy)
            m.Vz = np.dot(a, v2.Vz)
        else:
            m = np.dot(a, b)

        return m

    def rdivide(self, b):
        # RDIVIDE Right division(components)
        #
        # Vd = RDIVIDE(V, B) Right division(Vd=V. / B).
        # V is a Vector whose components Vx, Vy and Vz are the components
        # of V divided by the scalar( or scalar matrix) B.
        #
        # See also Vector.

        v_d = copy.deepcopy(self)
        v_d.Vx = np.divide(self.Vx, b)
        v_d.Vy = np.divide(self.Vy, b)
        v_d.Vz = np.divide(self.Vz, b)

        return v_d

    def versor(self):
        # VERSOR Unitary vector
        #
        # U = VERSOR(V) returns the unit vector set corresponding to
        # the vector set V.
        # The coordinates X, Y and Z of U are the ones of V.
        #
        # See also Vector.
        v = copy.deepcopy(self)
        v_temp = np.array([v.Vx, v.Vy, v.Vz])
        norm = np.linalg.norm(v_temp)
        if norm == 0:
            return v
        return v.rdivide(norm)

    def topoint(self):
        # TOPOINT Vector to point
        #
        # P = TOPOINT(V) converts the set of vectors V into the set of
        # points P.The coordinates X, Y and Z of the points are the
        # components of V.
        #
        # See also Vector, Point.
        v = copy.deepcopy(self)

        return Point(np.real(v.Vx), np.real(v.Vy), np.real(v.Vz))

    def toline(self):
        from SLine import SLine
        # TOLINE Vector to line
        #
        # LN = TOLINE(V) converts the set of vectors V into the set of
        # lines LN.The coordinates X, Y and Z of the initial points
        # of the lines are the coordinates of V and the coordinates
        # of the final points are the sum of the coordinates of V and
        # of the components of V.
        #
        # See also Vector, SLine.
        v = copy.deepcopy(self)

        return SLine(Point(v.X, v.Y, v.Z), Point(v.X + np.real(v.Vx), v.Y + np.real(v.Vy), v.Z + np.real(v.Vz)))
