import math
import numpy as np
import copy


class SLine:
    def __init__(self, p1, p2):
        '''
        Generates a line from two points
        :param p1: first point
        :param p2: second point
        '''
        self.p1 = p1
        self.p2 = p2

    def disp(self):
        '''
        Displays set of lines state
        '''
        print(self.p1.X)
        print(self.p1.Y)
        print(self.p1.Z)
        print(self.p2.X)
        print(self.p2.Y)
        print(self.p2.Z)

    def translate(self, dp):
        '''
        Translate set of lines, seld, by point or vector dp.
        If dP is a Point, the translation corresponds to the coordinates X, Y and Z.
        If dP is a Vector, the translation corresponds to the components Vx, Vy and Vz.
        :param dp: point of translation
        :return:
        '''
        ln = copy.deepcopy(self)
        ln_t = ln
        ln_t.p1 = ln_t.p1.translate(dp)
        ln_t.p2 = ln_t.p2.translate(dp)
        return ln_t

    def xrotation(self, phi):
        '''
        Rotates set of lines around x - axis
        :param phi: angle phi[rad] rotated around x - axis
        :return: new set of lines
        '''
        ln = copy.deepcopy(self)
        ln_r = ln
        ln_r.p1 = ln_r.p1.xrotation(phi)
        ln_r.p2 = ln_r.p2.xrotation(phi)
        return ln_r

    def yrotation(self, phi):
        '''
        Rotates set of lines around y - axis
        :param phi: angle phi[rad] rotated around y - axis
        :return: new set of lines
        '''
        ln = copy.deepcopy(self)
        ln_r = ln
        ln_r.p1 = ln_r.p1.yrotation(phi)
        ln_r.p2 = ln_r.p2.yrotation(phi)
        return ln_r

    def zrotation(self, phi):
        '''
        Rotates set of lines around z - axis
        :param phi: angle phi[rad] rotated around z - axis
        :return: new set of lines
        '''
        ln = copy.deepcopy(self)
        ln_r = ln
        ln_r.p1 = ln_r.p1.zrotation(phi)
        ln_r.p2 = ln_r.p2.zrotation(phi)
        return ln_r

    def numel(self):
        '''
        Calculates number of lines in set, self.
        :return: number of lines in set
        '''
        ln = copy.deepcopy(self)
        return ln.p1.size

    def size(self, varargin = None):
        '''
        Returns a two-element row vector with the number...
        of rows and columns in the line set LN.
        If argument is given returns the length of dimension...
        specified by the scalar DIM, input, in the line set, self.
        :param varargin:
        :return:
        '''
        ln = copy.deepcopy(self)
        if varargin is None:
            s = ln.p1.size()
        else:
            s = ln.p1.size(varargin[1])
        return s

    def angle(self,  ln2):
        '''
        Calculates the angle between the set of lines, self and ln2.
        :param ln2: set of lines
        :return: degrees between line
        '''
        ln1 = copy.deepcopy(self)
        phi = ln1.p2.minus(ln1.p1).angle(ln2.p2.minus(ln2.p1))
        return phi

    def versor(self):
        '''
        Returns unitvectors of lines in set
        :return: Returns unitvectors of lines in set
        '''
        from Vector import Vector
        ln = copy.deepcopy(self)
        v = Vector(np.zeros(ln.shape()), np.zeros(ln.shape()), np.zeros(ln.shape()), ln.p2.X - ln.p1.X,
                   ln.p2.Y - ln.p1.Y,
                   ln.p2.Z - ln.p1.Z)
        u = v.normalize()
        return u

    def tovector(self):
        '''
        Convert lines to vectors
        :return: set of vectors
        '''
        from Vector import Vector
        ln = copy.deepcopy(self)
        v = Vector(ln.p1.X, ln.p1.Y, ln.p1.Z, ln.p2.X - ln.p1.X, ln.p2.Y - ln.p1.Y, ln.p2.Z - ln.p1.Z)
        return v
