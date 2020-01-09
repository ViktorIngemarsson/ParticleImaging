from Point import Point
import numpy as np
import copy
from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.pyplot as plt
from SLine import SLine


class Vector(Point):
    def __init__(self, X, Y, Z, xcomponent, ycomponent, zcomponent):
        '''
        Generates a set of vectos
        :param X: Coordinates x
        :param Y: Coordinates y
        :param Z: Coordinates z
        :param xcomponent: Vector component x
        :param ycomponent: Vector component y
        :param zcomponent: Vector component z
        '''
        super().__init__(X, Y, Z)
        self.X = X
        self.Y = Y
        self.Z = Z
        self.Vx = xcomponent
        self.Vy = ycomponent
        self.Vz = zcomponent

    def disp(self):
        '''
        Display set of vectors
        '''
        print(self.X)
        print(self.Y)
        print(self.Z)
        print(self.Vx)
        print(self.Vy)
        print(self.Vz)

    def plot(self, varargin=None):
        '''
        Plot set of vectors
        '''
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.quiver(self.X, self.Y, self.Z, self.Vx, self.Vy, self.Vz, length=0.1, normalize=True)
        plt.show()

    def plot_multiple_vectors(self, moarVectors, radius):
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.quiver(self.X, self.Y, self.Z, self.Vx, self.Vy, self.Vz, length=radius, normalize=True)
        for vec in moarVectors:
            ax.quiver(vec.X, vec.Y, vec.Z, vec.Vx, vec.Vy, vec.Vz, length=radius, normalize=True)

        u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
        x = np.cos(u) * np.sin(v)
        y = np.sin(u) * np.sin(v)
        z = np.cos(v)
        ax.plot_surface(np.multiply(x, radius), np.multiply(y, radius), np.multiply(z, radius), rstride=1, cstride=1,
                       color="g", alpha=0.5, linewidth=0)
        plt.show()

        return

    def xrotation(self, phi):
        '''
        Rotates set of vectors V around x - axis
        :param phi: angle phi[rad] rotated around x-axis
        :return: new set of vectors
        '''
        v = copy.deepcopy(self)
        v_r = Point.xrotation(v, phi)
        v_r.Vy = v.Vy * np.cos(phi) - v.Vz * np.sin(phi)
        v_r.Vz = v.Vy * np.sin(phi) + v.Vz * np.cos(phi)
        return v_r

    def yrotation(self, phi):
        '''
        Rotates set of vectors V around y - axis
        :param phi: angle phi[rad] rotated around y-axis
        :return: new set of vectors
        '''
        v = copy.deepcopy(self)
        v_r = Point.yrotation(v, phi)
        v_r.Vx = v.Vx * np.cos(phi) + v.Vz * np.sin(phi)
        v_r.Vz = -v.Vx * np.sin(phi) + v.Vz * np.cos(phi)
        return v_r

    def zrotation(self, phi):
        '''
        Rotates set of vectors V around z - axis
        :param phi: angle phi[rad] rotated around z-axis
        :return: new set of vectors
        '''
        v = copy.deepcopy(self)
        v_r = Point.zrotation(v, phi)
        v_r.Vx = v.Vx * np.cos(phi) - v.Vy * np.sin(phi)
        v_r.Vy = v.Vx * np.sin(phi) + v.Vy * np.cos(phi)
        return v_r

    def uminus(self):
        '''
        Inverts the components Vx, Vy and Vz of V.
        :return: Inverted vector
        '''
        v_m = copy.deepcopy(self)
        v_m.Vx = -v_m.Vx
        v_m.Vy = -v_m.Vy
        v_m.Vz = -v_m.Vz

        return v_m

    def plus(self, v2):
        '''
        The components Vx, Vy and Vz of V are the sum of the components of V1 and V2.
        :param v2: vector to be added to self
        :return: new summarized vector
        '''
        v = copy.deepcopy(self)
        v.Vx = self.Vx + v2.Vx
        v.Vy = self.Vy + v2.Vy
        v.Vz = self.Vz + v2.Vz
        return v

    def minus(self, v2):
        '''
        The components Vx, Vy and Vz of V are the difference of the components of V1 and V2.
        :param v2: vector to be subtracted from self
        :return: new vector of difference in vector components
        '''
        v = copy.deepcopy(self)
        v.Vx = self.Vx - v2.Vx
        v.Vy = self.Vy - v2.Vy
        v.Vz = self.Vz - v2.Vz

        return v

    def mtimes(self, b):
        '''
        Computes vector product of self and b
        :param b: vector or scalar to compute product with self
        :return: vector/scalar product
        '''
        a = copy.deepcopy(self)
        if isinstance(a, Vector) & isinstance(b, Vector):
            v1 = copy.deepcopy(a)
            v2 = b
            m = copy.deepcopy(v1)
            m.Vx = np.multiply(v1.Vy, v2.Vz) - np.multiply(v1.Vz, v2.Vy)
            m.Vy = -np.multiply(v1.Vx, v2.Vz) + np.multiply(v1.Vz, v2.Vx)
            m.Vz = np.multiply(v1.Vx, v2.Vy) - np.multiply(v1.Vy, v2.Vx)
        elif isinstance(a, Vector):
            v1 = copy.deepcopy(a)
            m = copy.deepcopy(v1)
            m.X = m.X * np.ones([b.size()], dtype=int)
            m.Y = m.Y * np.ones([b.size()], dtype=int)
            m.Z = m.Z * np.ones([b.size()], dtype=int)
            m.Vx = np.multiply(v1.Vx, b)
            m.Vy = np.multiply(v1.Vy, b)
            m.Vz = np.multiply(v1.Vz, b)
        elif isinstance(b, Vector):
            v2 = copy.deepcopy(b)
            m = copy.deepcopy(v2)
            m.X = m.X * np.ones([a.size()], dtype=int)
            m.Y = m.Y * np.ones([a.size()], dtype=int)
            m.Z = m.Z * np.ones([a.size()], dtype=int)
            m.Vx = np.multiply(a, v2.Vx)
            m.Vy = np.multiply(a, v2.Vy)
            m.Vz = np.multiply(a, v2.Vz)
        else:
            m = np.multiply(a, b)

        return m

    def times(self, b):
        '''
        Computes scalar product of self and b
        :param b: vector or scalar
        :return:scalar matrix if b is vector, vector if b is scalar
        '''
        a = copy.deepcopy(self)

        if isinstance(a, Vector) & isinstance(b, Vector):
            v1 = copy.deepcopy(a)
            v2 = copy.deepcopy(b)
            m = np.multiply(v1.Vx, v2.Vx) + np.multiply(v1.Vy, v2.Vy) + np.multiply(v1.Vz, v2.Vz)
        elif isinstance(a, Vector):
            v1 = copy.deepcopy(a)
            m = copy.deepcopy(v1)
            m.X = m.X * np.ones([b.size()], dtype=int)
            m.Y = m.Y * np.ones([b.size()], dtype=int)
            m.Z = m.Z * np.ones([b.size()], dtype=int)
            m.Vx = np.multiply(v1.Vx, b)
            m.Vy = np.multiply(v1.Vy, b)
            m.Vz = np.multiply(v1.Vz, b)
        elif isinstance(b, Vector):
            v2 = copy.deepcopy(b)
            m = copy.deepcopy(v2)
            m.X = m.X * np.ones([a.size()], dtype=int)
            m.Y = m.Y * np.ones([a.size()], dtype=int)
            m.Z = m.Z * np.ones([a.size()], dtype=int)
            m.Vx = np.multiply(a, v2.Vx)
            m.Vy = np.multiply(a, v2.Vy)
            m.Vz = np.multiply(a, v2.Vz)
        else:
            m = np.multiply(a, b)

        return m

    def rdivide(self, b):
        '''
        Right division self/b
        :param b: scalar or scalar matrix
        :return: vector, self/b
        '''
        v_d = copy.deepcopy(self)
        v_d.Vx = np.divide(self.Vx, b)
        v_d.Vy = np.divide(self.Vy, b)
        v_d.Vz = np.divide(self.Vz, b)

        return v_d

    def versor(self):
        '''
        Return unit vector of self
        :return: unit vector of self
        '''
        v = copy.deepcopy(self)
        v_temp = np.transpose(np.array([v.Vx, v.Vy, v.Vz]))
        norm = np.linalg.norm(v_temp, axis=1)
        return v.rdivide(norm)

    def topoint(self):
        '''
        Converts set of vectors, self, to set of points
        :return: set of points
        '''
        v = copy.deepcopy(self)
        return Point(np.real(v.Vx), np.real(v.Vy), np.real(v.Vz))

    def toline(self):
        '''
        Converts set of vectors to set of lines.
        The coordinates X, Y and Z of the initial points ...
        of the lines are the coordinates of V and the...
        coordinates of the final points are the sum of the coordinates of V and
        :return: set of lines
        '''
        v = copy.deepcopy(self)
        return SLine(Point(v.X, v.Y, v.Z), Point(v.X + np.real(v.Vx), v.Y + np.real(v.Vy), v.Z + np.real(v.Vz)))
