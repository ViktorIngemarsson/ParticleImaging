import numpy as np
import matplotlib.pyplot as plt
from Vector import Vector
from SLine import SLine

class Point:
    def __init__(self, X, Y, Z):
        self.X = X
        self.Y = Y
        self.Z = Z

    def plot(self,varargin):
        # A 3d plot of the point (self), with the properties varargin.
        # TODO: Make the option varargin have an impact on the plotting
        ax = plt.axes(projection='3d')
        ax.scatter3D(self.X, self.Y, self.Z, c=self.Z, cmap='Greens');
        return

    def disp(self):
        print(self.X)
        print(self.Y)
        print(self.Z)
        return

    def translate(self,dp):
        p_t = self
        if  dp.isinstance(Vector):
            p_t.X = p_t.X + dp.Vx
            p_t.Y = p_t.Y + dp.Vy
            p_t.Z = p_t.Z + dp.Vz
        elif dp.isinstance(Point):
            p_t.X = p_t.X + dp.X
            p_t.Y = p_t.Y + dp.Y
            p_t.Z = p_t.Z + dp.Z
        return p_t

    def xrotation(self,phi):
        p_r = self
        p_r.Y = self.Y.dot(np.cos(phi)) - self.Z.dot(np.sin(phi))
        p_r.Z = self.Y.dot(np.sin(phi)) + self.Z.dot(np.cos(phi))
        return p_r

    def yrotation(self,phi):
        p_r = self
        p_r.X = self.X.dot(np.cos(phi)) + self.Z.dot(np.sin(phi))
        p_r.Z = -self.X.dot(np.sin(phi)) + self.Z.dot(np.cos(phi))
        return p_r

    def zrotation(self,phi):
        p_r = self
        p_r.X = self.X.dot(np.cos(phi)) - self.Y.dot(np.sin(phi))
        p_r.Y = self.X.dot(np.sin(phi)) + self.Y.dot(np.cos(phi))
        return p_r

    def numel(self):
        return np.size(self.X)

    def size(self,varargin):
        if not varargin:
            return np.shape(self.X, varargin[0])
        else:
            return np.shape(self.X)

    def uplus(self):
        return self

    def uminus(self):
        p_m = self
        p_m.X = -p_m.X
        p_m.Y = -p_m.Y
        p_m.Z = -p_m.Z
        return p_m

    def plus(self,p1,p2):
        p = p1
        p.X = p1.X + p2.X
        p.Y = p1.Y + p2.Y
        p.Z = p1.Z + p2.Z
        return p

    def minus(self,p1,p2):
        p = p1
        p.X = p1.X - p2.X
        p.Y = p1.Y - p2.Y
        p.Z = p1.Z - p2.Z
        return p

    def mtimes(self,a,b):
        if a.isinstance(Point) and b.isinstance(Point):
            p1 = a
            p2 = b
            m = p1
            m.X = p1.Y.dot(p2.Z) - p1.Z.dot(p2.Y)
            m.Y = -p1.X.dot(p2.Z) + p1.Z.dot(p2.X)
            m.Z = p1.X.dot(p2.Y) - p1.Y.dot(p2.X)
        elif a.isinstance(Point):
            p1 = a
            m = p1
            m.X = p1.X.dot(b)
            m.Y = p1.Y.dot(b)
            m.Z = p1.Z.dot(b)

        elif b.isinstance(Point):
            p2 = b
            m = p2
            m.X = a.dot(p2.X)
            m.Y = a.dot(p2.Y)
            m.Z = a.dot(p2.Z)
        else:
            m = a.dot(b)

        return m

    def times(self,a,b):
        if a.isinstance(Point) and b.isinstance(Point):
            p1 = a
            p2 = b
            m = p1.X.dot(p2.X) + p1.Y.dot(p2.Y) + p1.Z.dot(p2.Z)
        elif a.isinstance(Point):
            p1 = a
            m = p1
            m.X = p1.X.dot(b)
            m.Y = p1.Y.dot(b)
            m.Z = p1.Z.dot(b)

        elif b.isinstance(Point):
            p2 = b
            m = p2
            m.X = a.dot(p2.X)
            m.Y = a.dot(p2.Y)
            m.Z = a.dot(p2.Z)

        else:
            m = a.dot(b)
        return m

    def rdivide(self,p,b):
        p_d = p
        p_d.X = np.divide(p.X,b)
        p_d.Y = np.divide(p.Y,b)
        p_d.Z = np.divide(p.Z,b)
        return p_d

    def norm(self):
        return np.linalg.norm(self)

    def normalize(self):
        return np.divide(self,self.norm())

    def angle(self,p1,p2):
        p1 = p1.normalize()
        p2 = p2.normalize()
        return np.real(np.arccos(p1.dot(p2)))

    def toline(self):
        return SLine(Point(np.zeros(np.shape(self)),np.zeros(np.shape(self)),np.zeros(np.shape(self))),self)

    def tovector(self):
        return Vector(np.zeros(np.shape(self)),np.zeros(np.shape(self)),np.zeros(np.shape(self)),self.X,self.Y,self.Z)