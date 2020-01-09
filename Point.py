import numpy as np
import matplotlib.pyplot as plt
from SLine import SLine
import copy
from mpl_toolkits.mplot3d import axes3d, Axes3D


class Point:
    def __init__(self, X, Y, Z):
        '''
        Generates a set of points
        :param X: coordinates of set in x
        :param Y: coordinates of set in y
        :param Z: coordinates of set in z
        '''
        self.X = X
        self.Y = Y
        self.Z = Z

    def disp(self):
        '''
        Display state of the set of points
        '''
        print(self.X)
        print(self.Y)
        print(self.Z)

    def translate(self, dp):
        '''
        Translates set of points P by dP.
        If dP is a Point, the translation corresponds to the...
        coordinates X, Y and Z.
        If dP is a Vector, the translation corresponds to the...
        components Vx, Vy and Vz.
        If dP is neither a Point or a Vector returns an error.
        :param dp: point or vector
        :return: translated set of points
        '''
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
        '''
        Rotates set of points around x - axis
        :param phi: angle phi[rad] rotated around x - axis
        :return: new set of points
        '''
        p_r = copy.deepcopy(self)
        p_r.Y = self.Y * np.cos(phi) - self.Z * np.sin(phi)
        p_r.Z = self.Y * np.sin(phi) + self.Z * np.cos(phi)
        return p_r

    def yrotation(self, phi):
        '''
        Rotates set of points around y - axis
        :param phi: angle phi[rad] rotated around y - axis
        :return: new set of points
        '''
        p_r = copy.deepcopy(self)
        p_r.X = self.X * np.cos(phi) + self.Z * np.sin(phi)
        p_r.Z = -self.X * np.sin(phi) + self.Z * np.cos(phi)
        return p_r

    def zrotation(self, phi):
        '''
        Rotates set of points around z - axis
        :param phi: angle phi[rad] rotated around z - axis
        :return: new set of points
        '''
        p_r = copy.deepcopy(self)
        p_r.X = self.X * np.cos(phi) - self.Y * np.sin(phi)
        p_r.Y = self.X * np.sin(phi) + self.Y * np.cos(phi)
        return p_r

    def numel(self):
        '''
        Returns the number of points in set self.
        :return: Number of points in set
        '''
        return np.size(self.X)

    def size(self, varargin=None):
        '''
        Returns a two-element row vector with the number...
        of rows and columns in the point set P.
        If given argument returns the length of the dimension specified...
        by the scalar DIM in the point set P.
        :param varargin:
        :return:length of specified dimension
        '''
        if varargin != None:
            return np.shape(np.asarray(self.X), varargin[0])
        else:
            return np.shape(np.asarray(self.X))

    def uminus(self):
        '''
        Inverts all points in set self coordinates
        :return: inverted set of points
        '''
        p_m = copy.deepcopy(self)
        p_m.X = -p_m.X
        p_m.Y = -p_m.Y
        p_m.Z = -p_m.Z
        return p_m

    def plus(self, p2):
        '''
        Adds coordinates of set self to set p2 of points
        :param p2: set of points
        :return: new set of points
        '''
        p1 = copy.deepcopy(self)
        p1.X = p1.X + p2.X
        p1.Y = p1.Y + p2.Y
        p1.Z = p1.Z + p2.Z
        return p1

    def minus(self, p2):
        '''
        Subtracts coordinates of set self to set p2 of points
        :param p2: set of points
        :return: new set of points
        '''
        p1 = copy.deepcopy(self)
        p1.X = p1.X - p2.X
        p1.Y = p1.Y - p2.Y
        p1.Z = p1.Z - p2.Z
        return p1

    def mtimes(self, b):
        '''
        Vector product whose coordinates X, Y and Z are the vector
        product of the coordinates of self and b.
        If b is scalar then scalar product.
        :param b: a vector or a scalar to be mutiplied by
        :return: set of points
        '''
        a = copy.deepcopy(self)
        if isinstance(a, Point) and isinstance(b, Point):
            p1 = copy.deepcopy(a)
            p2 = b
            m = copy.deepcopy(p1)
            m.X = np.subtract(np.multiply(p1.Y, p2.Z), np.multiply(p1.Z, p2.Y))
            m.Y = np.add(np.multiply(-p1.X, p2.Z), np.multiply(p1.Z, p2.X))
            m.Z = np.subtract(np.multiply(p1.X, p2.Y), np.multiply(p1.Y, p2.X))
        elif isinstance(a, Point):
            p1 = copy.deepcopy(a)
            m = copy.deepcopy(p1)
            m.X = np.multiply(p1.X, b)
            m.Y = np.multiply(p1.Y, b)
            m.Z = np.multiply(p1.Z, b)
        elif isinstance(b, Point):
            p2 = copy.deepcopy(b)
            m = copy.deepcopy(p2)
            m.X = np.multiply(a, p2.X)
            m.Y = np.multiply(a, p2.Y)
            m.Z = np.multiply(a, p2.Z)
        else:
            m = np.multiply(a, b)

        return m

    def times(self, b):
        '''
        Calculates scalar product between set self and set of points b.
        If is scalar then regular scalar product.
        :param b: set of points or scalar
        :return: set of points
        '''
        a = copy.deepcopy(self)
        if isinstance(a, Point) and isinstance(b, Point):
            p1 = copy.deepcopy(a)
            p2 = b
            m = (p1.X * p2.X) + p1.Y * p2.Y + p1.Z * p2.Z #Changed but might have worked
        elif isinstance(a, Point):
            p1 = a
            m = p1
            m.X = np.multiply(p1.X, b)
            m.Y = np.multiply(p1.Y, b)
            m.Z = np.multiply(p1.Z, b)

        elif isinstance(b, Point):
            p2 = b
            m = p2
            m.X = np.multiply(a, p2.X)
            m.Y = np.multiply(a, p2.Y)
            m.Z = np.multiply(a, p2.Z)

        else:
            m = np.multiply(a, b)
        return m

    def rdivide(self, b):
        '''
        Right scalar division between set self and scalar (or scalar matrix) b.
        :param b: scalar or scalar matrix
        :return: set of points
        '''
        p_d = copy.deepcopy(self)
        p_d.X = np.divide(self.X, b)
        p_d.Y = np.divide(self.Y, b)
        p_d.Z = np.divide(self.Z, b)
        return p_d

    def norm(self):
        '''
        Calculates norm of set of points
        :return: norm of set of points
        '''
        p = copy.deepcopy(self)
        return np.sqrt(p.times(p))

    def normalize(self):
        '''
        Normalize the set of points
        :return: set of points
        '''
        nominator = copy.deepcopy(self)
        return nominator.rdivide(nominator.norm())

    def angle(self, p2):
        '''
        Calculates angle between set of points self and set p2
        :param p2: set of points
        :return: angle
        '''
        p1 = copy.deepcopy(self)
        p1 = p1.normalize()
        p2 = p2.normalize()
        return np.real(np.arccos(p1.times(p2)))

    def toline(self):
        '''
        Convert set of points to srt of lines
        :return: set of lines
        '''
        temp = copy.deepcopy(self)
        return SLine(Point(np.zeros(temp.size()), np.zeros(temp.size()), np.zeros(temp.size())), temp)

    def tovector(self):
        '''
        Convert set of points to set of vectors
        :return: set of vectors
        '''
        from Vector import Vector
        temp = copy.deepcopy(self)
        return Vector(np.zeros(temp.size()), np.zeros(temp.size()), np.zeros(temp.size()), temp.X, temp.Y,
                      temp.Z)

