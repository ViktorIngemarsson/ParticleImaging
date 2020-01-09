import numpy as np
import matplotlib.pyplot as plt
from Vector import Vector
from Spherical import Spherical
from Point import Point
from Plane import Plane
from Transform import Transform
from operator import itemgetter
import math
import copy


class Ray:
    def __init__(self, v, P, pol):
        '''
        Generates set of rays
        :param v: vectors
        :param P: power of rays
        :param pol: polarisation of rays
        '''
        self.v = v
        self.P = P
        self.pol = pol
        self.pol.X = v.X
        self.pol.Y = v.Y
        self.pol.Z = v.Z

    def disp(self):
        '''
        Displays set of ray state state
        '''
        print(self.v.X)
        print(self.v.Y)
        print(self.v.Z)
        print(self.v.Vx)
        print(self.v.Vy)
        print(self.v.Vz)
        print(self.P)
        print(self.pol.Vx)
        print(self.pol.Vy)
        print(self.pol.Vz)

    def translate(self, dp):
        '''
        Translates set of rays R by dP.
        If dP is a Point, the translation corresponds to the...
        coordinates X, Y and Z.
        If dP is a Vector, the translation corresponds to the...
        components Vx, Vy and Vz.
        :param dp:
        :return:
        '''
        r = copy.deepcopy(self)
        r_t = r
        r_t.v = r.v.translate(dp)
        r_t.pol = r.pol.translate(dp)
        return r_t

    def xrotation(self, phi):
        '''
        Rotates set of rays around x - axis
        :param phi: angle phi[rad] rotated around x-axis
        :return: new set of vectors
        '''
        r = copy.deepcopy(self)
        r_r = r
        r_r.v = r.v.xrotation(phi)
        r_r.pol = r.pol.xrotation(phi)
        return r_r

    def yrotation(self, phi):
        '''
        Rotates set of rays around y - axis
        :param phi: angle phi[rad] rotated around y - axis
        :return: new set of vectors
        '''
        r = copy.deepcopy(self)
        r_r = r
        r_r.v = r.v.yrotation(phi)
        r_r.pol = r.pol.yrotation(phi)
        return r_r

    def zrotation(self, phi):
        '''
        Rotates set of rays around z - axis
        :param phi: angle phi[rad] rotated around z - axis
        :return: new set of vectors
        '''
        r = copy.deepcopy(self)
        r_r = r
        r_r.v = r.v.zrotation(phi)
        r_r.pol = r.pol.zrotation(phi)
        return r_r

    def numel(self):
        '''
        Gives size of set of rays
        :return: number of rays in set
        '''
        r = copy.deepcopy(self)
        return np.size(r.v.X)

    def size(self, varargin=None):
        '''
        Returns a two-element row vector with the number ...
        of rows and columns in the ray set R.
        If given argument returns the length of the dimension specified
        by the scalar DIM in the ray set self.
        :param varargin:
        :return: size of ray set
        '''
        if varargin != None:
            return np.shape(np.asarray(self.v.X), varargin[0])
        else:
            return np.shape(np.asarray(self.v.X))

    def uminus(self):
        '''
        Returns ray with inverted vector values
        :return: inverted set of rays
        '''
        r_m = copy.deepcopy(self)
        r_m.v = r_m.v.uminus()
        return r_m

    def angle(self,r2):
        '''
        Calculates angle between set of ray, self and r2
        :param r2: set of rays
        :return: angle between set of rays
        '''
        r1 = copy.deepcopy(self)
        return r1.v.angle(r2.v)

    def versor(self):
        '''
        Calculates unitvector of set of rays
        :return: Set of unit vectors
        '''
        return copy.deepcopy(self).v.versor()

    def toline(self):
        '''
        Convert set of rays to set of lines
        :return: set of lines
        '''
        return copy.deepcopy(self).v.toline()

    def snellslaw(self, s, n1, n2, n):
        '''
        Calculates the reflection and refraction between the set of rays, self,...
        and the spherical object s. Where the rays comes from medium 1 with...
        refractive index n1 into spherical object with refractive index...
        of n2. Calculates what happens looking at the n:th intersection between the two.
        :param s: Particle spherical object
        :param n1: refractive index that rays are coming from
        :param n2: refractive index of the spherical object
        :param n: what index of intersection to look at
        :return:
        '''
        r = copy.deepcopy(self)
        if isinstance(s, Plane):
            pl = s
            # intersection between ray and plane
            p = pl.intersectionpoint(r)

            # normal line
            perp = pl.perpline(p)

            # Incidence angle
            theta_i = r.toline().angle(perp)
            if (theta_i > np.pi / 2).any():
                theta_i[theta_i > np.pi / 2] = -theta_i[theta_i > np.pi / 2] + np.pi

            # Transmission angle
            theta_t = np.arcsin(n1 / n2 * np.sin(theta_i))

            # translation to origin(0)
            r0 = r.translate(p.uminus())
            perp0 = perp.translate(p.uminus())

            # rotation around z(1)
            phi1 = np.arctan2(perp0.p2.X, perp0.p2.Y)
            r1 = r0.zrotation(phi1)
            perp1 = perp0.zrotation(phi1)

            # rotation around x(2)
            phi2 = np.arctan2(perp1.p2.Y, perp1.p2.Z)
            r2 = r1.xrotation(phi2)

            # rotation around z(3)
            phi3 = np.arctan2(r2.v.X, r2.v.Y)
            r3 = r2.zrotation(phi3)

            output_dictS = Ray.rtcoeffs(theta_i, n1, n2)
            Rs = output_dictS["Rs"]
            Ts = output_dictS["Ts"]
            rs = output_dictS["rs"]
            ts = output_dictS["ts"]

            output_dictP = Ray.rtcoeffp(theta_i, n1, n2)
            Rp = output_dictP["Rp"]
            Tp = output_dictP["Tp"]
            rp = output_dictP["rp"]
            tp = output_dictP["tp"]

            dtheta = theta_i *2
            dtheta[r3.v.Z < 0] = -dtheta[r3.v.Z < 0]
            r_r3 = r3.xrotation(dtheta).uminus()

            r_r3.pol.Vx = np.abs(rs * r_r3.pol.Vx)
            r_r3.pol.Vy = np.abs(rp * r_r3.pol.Vy)
            r_r3.pol.Vz = np.abs(rp * r_r3.pol.Vz)

            cs = np.divide((np.square(r3.pol.Vx)), np.square(r3.pol.norm()))
            cp = np.divide((np.square(r3.pol.Vy)) + np.square(r3.pol.Vz), np.square(r3.pol.norm()))
            R = cs * Rs + cp * Rp

            r_r3.v.X = np.zeros(np.shape(r))
            r_r3.v.Y = np.zeros(np.shape(r))
            r_r3.v.Z = np.zeros(np.shape(r))
            r_r3.P = r.P * R
            r_r3.pol.X = np.zeros(np.shape(r))
            r_r3.pol.Y = np.zeros(np.shape(r))
            r_r3.pol.Z = np.zeros(np.shape(r))

            dtheta = -np.pi + theta_i - theta_t
            dtheta[r3.v.Z < 0] = -dtheta[r3.v.Z < 0]
            r_t3 = r3.xrotation(np.real(dtheta)).uminus()
            r_t3.pol.Vy = -r_t3.pol.Vy
            r_t3.pol.Vz = -r_t3.pol.Vz

            r_t3.pol.Vx = ts * r_t3.pol.Vx
            r_t3.pol.Vy = tp * r_t3.pol.Vy
            r_t3.pol.Vz = tp * r_t3.pol.Vz

            r_t3.pol.Vx = np.multiply(ts, r_t3.pol.Vx)
            r_t3.pol.Vy = np.multiply(tp, r_t3.pol.Vy)
            r_t3.pol.Vz = np.multiply(tp, r_t3.pol.Vz)

            T = 1 - R

            r_t3.v.X = np.zeros(np.shape(r))
            r_t3.v.Y = np.zeros(np.shape(r))
            r_t3.v.Z = np.zeros(np.shape(r))
            r_t3.P = np.multiply(r.P, T)
            r_t3.pol.X = np.zeros(np.shape(r))
            r_t3.pol.Y = np.zeros(np.shape(r))
            r_t3.pol.Z = np.zeros(np.shape(r))

            tir = np.imag(theta_t) != 0

            r_r3.P[tir] = r.P[tir]
            r_r3.pol.Vx[tir] = np.abs(r_r3.pol.Vx[tir])
            r_r3.pol.Vy[tir] = np.abs(r_r3.pol.Vy[tir])
            r_r3.pol.Vz[tir] = np.abs(r_r3.pol.Vz[tir])

            r_t3.P[tir] = 0
            r_t3.v.Vx[tir] = np.real(r_t3.v.Vx[tir])
            r_t3.v.Vy[tir] = np.real(r_t3.v.Vy[tir])
            r_t3.v.Vz[tir] = np.real(r_t3.v.Vz[tir])
            r_t3.pol.Vx[tir] = np.abs(r_t3.pol.Vx[tir])
            r_t3.pol.Vy[tir] = np.abs(r_t3.pol.Vy[tir])
            r_t3.pol.Vz[tir] = np.abs(r_t3.pol.Vz[tir])

            # back rotation around z(-3)
            r_r4 = r_r3.zrotation(-phi3)
            r_t4 = r_t3.zrotation(-phi3)

            # back rotation around x(-2)
            r_r5 = r_r4.xrotation(-phi2)
            r_t5 = r_t4.xrotation(-phi2)

            # back rotation around z(-1)
            r_r6 = r_r5.zrotation(-phi1)
            r_t6 = r_t5.zrotation(-phi1)

            # back translation from origin (-0)
            r_r = r_r6.translate(p)
            r_t = r_t6.translate(p)

        elif isinstance(s, Spherical):
            # intersection between ray and sphere
            p = s.intersectionpoint(r, n)

            # normal line
            ln = s.perpline(p)

            # tangent plane
            pl = Plane.perpto(ln, p)

            # Snell's law
            return r.snellslaw(pl, n1, n2, 1)

        return r_r, r_t, perp

    @staticmethod
    def rtcoeffs(theta_i, n1, n2):
        '''
        calculates the Fresnel coefficient for an s-polarized ray impinging with...
        incidence angle theta_i on a planar surface. The refractive...
        indices of the incoming medium is n1 and the one of the...
        transmission medium n2.
        :param theta_i:
        :param n1:
        :param n2:
        :return Rs: Fresnel reflection coefficient for power
        :return Ts: Fresnel transmission coefficient for power
        :return rs: Fresnel reflection coefficient for electric field
        :return ts: Fresnel transmission coefficient for electric field
        '''
        theta_t = np.arcsin(n1 / n2 * np.sin(theta_i))

        rs = np.divide((n1 * np.cos(theta_i) - n2 * np.cos(theta_t)), (n1 * np.cos(theta_i) + n2 * np.cos(theta_t)))
        ts = np.divide(2 * n1 * np.cos(theta_i), (n1 * np.cos(theta_i) + n2 * np.cos(theta_t)))

        Rs = np.abs(np.square(rs))
        Ts = 1 - Rs

        return {
            "Rs": Rs,
            "Ts": Ts,
            "rs": rs,
            "ts": ts
        }

    @staticmethod
    def rtcoeffp(theta_i, n1, n2):
        '''
        calculates the Fresnel coefficient for a p-polarized ray impinging with
        incidence angle theta_i on a planar surface. The refractive
        indices of the incoming medium is n1 and the one of the
        transmission medium n2.
        :param theta_i:
        :param n1:
        :param n2:
        :return Rp: Fresnel reflection coefficient for power
        :return Tp: Fresnel transmission coefficient for power
        :return rp: Fresnel reflection coefficient for electric field
        :return tp: Fresnel transmission coefficient for electric field
        '''
        theta_t = np.arcsin(n1 / n2 * np.sin(theta_i))

        rp = np.divide((n1 * np.cos(theta_t) - n2 * np.cos(theta_i)), (n1 * np.cos(theta_t) + n2 * np.cos(theta_i)))
        tp = np.divide(2 * n1 * np.cos(theta_i), (n1 * np.cos(theta_t) + n2 * np.cos(theta_i)))

        Rp = np.abs(np.square(rp))
        Tp = 1 - Rp

        return {
            "Rp": Rp,
            "Tp": Tp,
            "rp": rp,
            "tp": tp
        }
