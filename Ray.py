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
        # v       directions (Vector)
        # P       powers (matrix)
        # pol     polarizations (Vector)
        self.v = v
        self.P = P
        self.pol = pol
        self.pol.X = v.X
        self.pol.Y = v.Y
        self.pol.Z = v.Z

    def disp(self):
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
        r = copy.deepcopy(self)
        r_t = r
        r_t.v = r.v.translate(dp)
        r_t.pol = r.pol.translate(dp)
        return r_t

    def xrotation(self, phi):
        r = copy.deepcopy(self)
        r_r = r
        r_r.v = r.v.xrotation(phi)
        r_r.pol = r.pol.xrotation(phi)
        return r_r

    def yrotation(self, phi):
        r = copy.deepcopy(self)
        r_r = r
        r_r.v = r.v.yrotation(phi)
        r_r.pol = r.pol.yrotation(phi)
        return r_r

    def zrotation(self, phi):
        r = copy.deepcopy(self)
        r_r = r
        r_r.v = r.v.zrotation(phi)
        r_r.pol = r.pol.zrotation(phi)
        return r_r

    def numel(self):
        r = copy.deepcopy(self)
        return np.size(r.v.X)

    def size(self, varargin=None):
        if varargin != None:
            return np.shape(np.asarray(self.v.X), varargin[0])
        else:
            return np.shape(np.asarray(self.v.X))

    def uplus(self):
        return copy.deepcopy(self)

    def uminus(self):
        r_m = copy.deepcopy(self)
        r_m.v = r_m.v.uminus()
        return r_m

    def angle(self,r2):
        r1 = copy.deepcopy(self)
        return r1.v.angle(r2.v)

    def versor(self):
        return copy.deepcopy(self).v.versor()

    def toline(self):
        return copy.deepcopy(self).v.toline()


    # TODO: I matlab så används nargin för att avgöra hur många input-argument functionen tar emot, hur gör vi det i python,
    #   och kan de olika funktionerna ta emot olika många inputs?

    def snellslaw(self, s, n1, n2, n):
        r = copy.deepcopy(self)
        if isinstance(s, Plane):
            pl = s
            # intersection between ray and plane
            p = pl.intersectionpoint(r)

            # normal line
            perp = pl.perpline(p)

            # Incidence angle
            theta_i = r.toline().angle(perp)
            if theta_i > np.pi / 2:
                theta_i = np.pi - theta_i

            # Transmission angle
            theta_t = np.arcsin(n1 / n2 * np.sin(theta_i))

            # translation to origin(0)
            r0 = r.translate(p.uminus())
            # pl0 = pl.translate(-p); %
            perp0 = perp.translate(p.uminus())
            # p0 = p.translate(-p); %

            # rotation around z(1)
            phi1 = np.arctan2(perp0.p2.X, perp0.p2.Y)
            r1 = r0.zrotation(phi1)
            # pl1 = pl0.zrotation(phi1); %
            perp1 = perp0.zrotation(phi1)
            # p1 = p0.zrotation(phi1); %

            # rotation around x(2)
            phi2 = np.arctan2(perp1.p2.Y, perp1.p2.Z)
            r2 = r1.xrotation(phi2)
            # pl2 = pl1.xrotation(phi2); %
            # perp2 = perp1.xrotation(phi2); %
            # p2 = p1.xrotation(phi2); %

            # rotation around z(3)
            phi3 = np.arctan2(r2.v.X, r2.v.Y)
            r3 = r2.zrotation(phi3)
            # pl3 = pl2.zrotation(phi3); %
            # perp3 = perp2.zrotation(phi3); %
            # p3 = p2.zrotation(phi3); %

            # Rs and Ts(appliesto r3.pol.Vx)
            output_dictS = Ray.rtcoeffs(theta_i, n1, n2)
            Rs = output_dictS["Rs"]
            Ts = output_dictS["Ts"]
            rs = output_dictS["rs"]
            ts = output_dictS["ts"]

            # Rp and Tp(r3.pol.Vy and Vz)
            output_dictP = Ray.rtcoeffp(theta_i, n1, n2)
            Rp = output_dictP["Rp"]
            Tp = output_dictP["Tp"]
            rp = output_dictP["rp"]
            tp = output_dictP["tp"]

            # Reflected ray
            if r3.v.Z > 0:
                r_r3 = r3.xrotation(2 * theta_i).uminus()
            else:
                r_r3 = r3.xrotation(-2 * theta_i).uminus()
            # dtheta = theta_i *2
            # if dtheta < 0: dtheta = -dtheta
            # r_r3 = r3.xrotation(dtheta).uminus()

            r_r3.pol.Vx = np.abs(np.dot(rs, r_r3.pol.Vx))
            r_r3.pol.Vy = np.abs(np.dot(rp, r_r3.pol.Vy))
            r_r3.pol.Vz = np.abs(np.dot(rp, r_r3.pol.Vz))

            cs = np.divide((np.square(r3.pol.Vx)), np.square(r3.pol.norm()))
            cp = np.divide((np.square(r3.pol.Vy)) + np.square(r3.pol.Vz), np.square(r3.pol.norm()))
            R = np.dot(cs, Rs) + np.dot(cp, Rp)

            r_r3.v.X = np.zeros(np.shape(r))
            r_r3.v.Y = np.zeros(np.shape(r))
            r_r3.v.Z = np.zeros(np.shape(r))
            r_r3.P = np.dot(r.P, R)
            r_r3.pol.X = np.zeros(np.shape(r))
            r_r3.pol.Y = np.zeros(np.shape(r))
            r_r3.pol.Z = np.zeros(np.shape(r))

            # transmitted ray
            if r3.v.Z > 0:
                r_t3 = r3.xrotation(-np.pi + theta_i - theta_t).uminus()
                r_t3.pol.Vy = -r_t3.pol.Vy
                r_t3.pol.Vz = -r_t3.pol.Vz
            else:
                r_t3 = r3.xrotation(np.pi - theta_i + theta_t).uminus()
                r_t3.pol.Vy = -r_t3.pol.Vy
                r_t3.pol.Vz = -r_t3.pol.Vz

            r_t3.pol.Vx = np.dot(ts, r_t3.pol.Vx)
            r_t3.pol.Vy = np.dot(tp, r_t3.pol.Vy)
            r_t3.pol.Vz = np.dot(tp, r_t3.pol.Vz)

            T = 1 - R

            r_t3.v.X = np.zeros(np.shape(r))
            r_t3.v.Y = np.zeros(np.shape(r))
            r_t3.v.Z = np.zeros(np.shape(r))
            r_t3.P = np.dot(r.P, T)
            r_t3.pol.X = np.zeros(np.shape(r))
            r_t3.pol.Y = np.zeros(np.shape(r))
            r_t3.pol.Z = np.zeros(np.shape(r))

            # ManagesTIR
            # indices = [1, 2, 3]
            # tir = [j for (i, j) in zip(np.imag(theta_t), indices) if i != 0]
            # tir = [x - 1 for x in tir]
            if np.imag(theta_t) != 0:
                # tir = np.imag(theta_t) != 0
                r_r3.P = r.P
                r_r3.pol.Vx = np.abs(r_r3.pol.Vx)
                r_r3.pol.Vy = np.abs(r_r3.pol.Vy)
                r_r3.pol.Vz = np.abs(r_r3.pol.Vz)

                r_t3.P = 0
                r_t3.v.Vx = np.real(r_t3.v.Vx)
                r_t3.v.Vy = np.real(r_t3.v.Vy)
                r_t3.v.Vz = np.real(r_t3.v.Vz)
                r_t3.pol.Vx = np.abs(r_t3.pol.Vx)
                r_t3.pol.Vy = np.abs(r_t3.pol.Vy)
                r_t3.pol.Vz = np.abs(r_t3.pol.Vz)

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
            output = r.snellslaw(pl, n1, n2, 1)
            r_r = output['r_r']
            r_t = output['r_t']
            perp = output['perp']

        return {
            "r_r": r_r,
            "r_t": r_t,
            "perp": perp
        }

    @staticmethod
    def rtcoeffs(theta_i, n1, n2):
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

    @staticmethod
    def beam2rays(b):

        X = np.dot(b.r, np.cos(b.phi))
        Y = b.r.dot(np.sin(b.phi))
        Z = np.zeros(np.shape(b.r))
        Vx = np.zeros(np.shape(b.r))
        Vy = np.zeros(np.shape(b.r))
        Vz = np.ones(np.shape(b.r))
        v = Vector(np.array(X, X), np.array(Y, Y), np.array(Z, Z), np.array(Vx, Vx), np.array(Vy, Vy), np.array(Vz, Vz))

        output_dict = Transform.Pol2CarVector(b.phi, b.Ephi, b.Er)
        Ex = output_dict["Ex"]
        Ey = output_dict["Ey"]

        pol = Vector(np.array(X, X), np.array(Y, Y), np.array(Z, Z), np.array(real(Ex), imag(Ex)),
                     np.array(real(Ey), imag(Ey)), np.array(np.zeros(np.shape(b.r)), np.zeros(np.shape(b.r))))

        dr = b.r(1, 2) - b.r(1, 1)
        dphi = b.phi(2, 1) - b.phi(1, 1)
        P = np.dot(b.intensity(), b.r) * dr * dphi
        P = np.array(np.divide(np.dot(P, np.square(np.abs(Ex))), (np.square(np.abs(Ex)) + np.square(np.abs(Ey)))),
                     np.divide(np.dot(P, np.square(np.abs(Ey))), (np.square(np.abs(Ex)) + np.square(np.abs(Ey)))))

        res = Ray(v, P, pol)

        return res

    @staticmethod
    def beam2focus(b, f):
        res = Ray.beam2rays(b)

        r = np.sqrt(np.square(res.v.X) + np.square(res.v.Y))
        theta = np.arcsin(r / f)

        res = res.translate(Point(0 * r, 0 * r, -f * np.cos(theta)))

        phi = np.arctan2(res.v.Y, res.v.X)
        res = res.zrotation(-phi)

        dp = Point(res.v.X, res.v.Y, res.v.Z)
        res = res.translate(-dp)

        res = res.yrotation(-theta)

        res = res.translate(dp)

        res = res.zrotation(phi)

        return res
