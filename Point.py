import numpy as np
import matplotlib.pyplot as plt
from Vector import Vector
from SLine import SLine
import copy


class Point:
    def __init__(self, X, Y, Z):
        self.X = X
        self.Y = Y
        self.Z = Z

    def plot(self, varargin):
        # A 3d plot of the point (self), with the properties varargin.
        # TODO: Make the option varargin have an impact on the plotting
        ax = plt.axes(projection='3d')
        ax.scatter3D(self.X, self.Y, self.Z, c=self.Z, cmap='Greens');
        return

    def disp(self):
        print(self.X)
        print(self.Y)
        print(self.Z)

    def translate(self, dp):
        p_t = copy.deepcopy(self)
        if isinstance(dp, Vector):
            p_t.X = p_t.X + dp.Vx
            p_t.Y = p_t.Y + dp.Vy
            p_t.Z = p_t.Z + dp.Vz
        elif isinstance(dp, Point):
            p_t.X = p_t.X + dp.X
            p_t.Y = p_t.Y + dp.Y
            p_t.Z = p_t.Z + dp.Z
        return p_t

    def xrotation(self, phi):
        p_r = copy.deepcopy(self)
        p_r.Y = np.dot(self.Y, np.cos(phi)) - np.dot(self.Z, np.sin(phi))
        p_r.Z = np.dot(self.Y, np.sin(phi)) + np.dot(self.Z, np.cos(phi))
        return p_r

    def yrotation(self, phi):
        p_r = copy.deepcopy(self)
        p_r.X = np.dot(self.X, np.cos(phi)) + np.dot(self.Z, np.sin(phi))
        p_r.Z = np.dot(-self.X, np.sin(phi)) + np.dot(self.Z, np.cos(phi))
        return p_r

    def zrotation(self, phi):
        p_r = copy.deepcopy(self)
        p_r.X = np.dot(self.X, np.cos(phi)) - np.dot(self.Y, np.sin(phi))
        p_r.Y = np.dot(self.X, np.sin(phi)) + np.dot(self.Y, np.cos(phi))
        return p_r

    def numel(self):
        return np.size(self.X)

    def size(self, varargin=None):
        if varargin != None:
            return np.size(self.X, varargin[0])

        else:
            return np.size(self.X)

    def uplus(self):
        return self

    def uminus(self):
        p_m = copy.deepcopy(self)
        p_m.X = -p_m.X
        p_m.Y = -p_m.Y
        p_m.Z = -p_m.Z
        return p_m

    def plus(self, p1, p2):
        p = p1
        p.X = p1.X + p2.X
        p.Y = p1.Y + p2.Y
        p.Z = p1.Z + p2.Z
        return p

    def minus(self, p1, p2):
        p = p1
        p.X = p1.X - p2.X
        p.Y = p1.Y - p2.Y
        p.Z = p1.Z - p2.Z
        return p

    def mtimes(self, a, b):
        if a.isinstance(Point) and b.isinstance(Point):
            p1 = a
            p2 = b
            m = p1
            m.X = np.dot(p1.Y, p2.Z) - np.dot(p1.Z, p2.Y)
            m.Y = np.dot(-p1.X, p2.Z) + np.dot(p1.Z, p2.X)
            m.Z = np.dot(p1.X, p2.Y) - np.dot(p1.Y, p2.X)
        elif a.isinstance(Point):
            p1 = a
            m = p1
            m.X = np.dot(p1.X, b)
            m.Y = np.dot(p1.Y, b)
            m.Z = np.dot(p1.Z, b)

        elif b.isinstance(Point):
            p2 = b
            m = p2
            m.X = np.dot(a, p2.X)
            m.Y = np.dot(a, p2.Y)
            m.Z = np.dot(a, p2.Z)
        else:
            m = np.dot(a, b)

        return m

    def times(self, a, b):
        if a.isinstance(Point) and b.isinstance(Point):
            p1 = a
            p2 = b
            m = np.dot(p1.X, p2.X) + np.dot(p1.Y, p2.Y) + np.dot(p1.Z, p2.Z)
        elif a.isinstance(Point):
            p1 = a
            m = p1
            m.X = np.dot(p1.X, b)
            m.Y = np.dot(p1.Y, b)
            m.Z = np.dot(p1.Z, b)

        elif b.isinstance(Point):
            p2 = b
            m = p2
            m.X = np.dot(a, p2.X)
            m.Y = np.dot(a, p2.Y)
            m.Z = np.dot(a, p2.Z)

        else:
            m = np.dot(a, b)
        return m

    def rdivide(self, p, b):
        p_d = p
        p_d.X = np.divide(p.X, b)
        p_d.Y = np.divide(p.Y, b)
        p_d.Z = np.divide(p.Z, b)
        return p_d

    def norm(self):
        return np.linalg.norm(self)

    def normalize(self):
        return np.divide(self, self.norm())

    def angle(self, p1, p2):
        p1 = p1.normalize()
        p2 = p2.normalize()
        return np.real(np.arccos(p1.dot(p2)))

    def toline(self):
        return SLine(Point(np.zeros(np.shape(self)), np.zeros(np.shape(self)), np.zeros(np.shape(self))), self)

    def tovector(self):
        return Vector(np.zeros(np.shape(self)), np.zeros(np.shape(self)), np.zeros(np.shape(self)), self.X, self.Y,
                      self.Z)
