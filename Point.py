import numpy as np
import matplotlib.pyplot as plt
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
        from Vector import Vector
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
            return np.shape(np.asarray(self.X), varargin[0])
        else:
            return np.shape(np.asarray(self.X))

    def uplus(self):
        return self

    def uminus(self):
        p_m = copy.deepcopy(self)
        p_m.X = -p_m.X
        p_m.Y = -p_m.Y
        p_m.Z = -p_m.Z
        return p_m

    def plus(self, p2):
        p1 = copy.deepcopy(self)
        p1.X = p1.X + p2.X
        p1.Y = p1.Y + p2.Y
        p1.Z = p1.Z + p2.Z
        return p1

    def minus(self, p2):
        p1 = copy.deepcopy(self)
        p1.X = p1.X - p2.X
        p1.Y = p1.Y - p2.Y
        p1.Z = p1.Z - p2.Z
        return p1

    def mtimes(self, b):
        a = copy.deepcopy(self)
        if isinstance(a, Point) and isinstance(b, Point):
            p1 = copy.deepcopy(a)
            p2 = b
            m = copy.deepcopy(p1)
            m.X = np.subtract(np.dot(p1.Y, p2.Z), np.dot(p1.Z, p2.Y))
            m.Y = np.add(np.dot(-p1.X, p2.Z), np.dot(p1.Z, p2.X))
            m.Z = np.subtract(np.dot(p1.X, p2.Y), np.dot(p1.Y, p2.X))
        elif isinstance(a, Point):
            p1 = copy.deepcopy(a)
            m = copy.deepcopy(p1)
            m.X = np.dot(p1.X, b)
            m.Y = np.dot(p1.Y, b)
            m.Z = np.dot(p1.Z, b)
        elif isinstance(b, Point):
            p2 = copy.deepcopy(b)
            m = copy.deepcopy(p2)
            m.X = np.dot(a, p2.X)
            m.Y = np.dot(a, p2.Y)
            m.Z = np.dot(a, p2.Z)
        else:
            m = np.dot(a, b)

        return m

    def times(self, b):
        a = copy.deepcopy(self)
        if isinstance(a, Point) and isinstance(b, Point):
            p1 = copy.deepcopy(a)
            p2 = b
            m = np.dot(p1.X, p2.X) + np.dot(p1.Y, p2.Y) + np.dot(p1.Z, p2.Z)
        elif isinstance(a, Point):
            p1 = a
            m = p1
            m.X = np.dot(p1.X, b)
            m.Y = np.dot(p1.Y, b)
            m.Z = np.dot(p1.Z, b)

        elif isinstance(b, Point):
            p2 = b
            m = p2
            m.X = np.dot(a, p2.X)
            m.Y = np.dot(a, p2.Y)
            m.Z = np.dot(a, p2.Z)

        else:
            m = np.dot(a, b)
        return m

    def rdivide(self, b):
        p_d = copy.deepcopy(self)
        p_d.X = np.divide(self.X, b)
        p_d.Y = np.divide(self.Y, b)
        p_d.Z = np.divide(self.Z, b)
        return p_d

    def norm(self):
        return np.linalg.norm(np.array([self.X, self.Y, self.Z]))

    def normalize(self):
        nominator = copy.deepcopy(self)
        return nominator.rdivide(nominator.norm())

    def angle(self, p2):
        p1 = copy.deepcopy(self)
        p1 = p1.normalize()
        p2 = p2.normalize()
        return np.real(np.arccos(p1.times(p2)))

    def toline(self):
        temp = copy.deepcopy(self)
        return SLine(Point(np.zeros(temp.size()), np.zeros(temp.size()), np.zeros(temp.size())), temp)

    def tovector(self):
        from Vector import Vector
        temp = copy.deepcopy(self)
        return Vector(np.zeros(temp.size()), np.zeros(temp.size()), np.zeros(temp.size()), temp.X, temp.Y,
                      temp.Z)
