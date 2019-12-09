import numpy as np
import matplotlib.pyplot as plt
from Vector import Vector
from SLine import SLine
from Point import Point
from Plane import Plane
from Transform import Transform
from operator import itemgetter
import math

def dotproduct(v1, v2):
    return sum((a * b) for a, b in zip(v1, v2))

def length(v):
    return math.sqrt(dotproduct(v, v))


def angle(v1, v2):
    return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

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

    # TODO: I matlab så används nargin för att avgöra hur många input-argument functionen tar emot, hur gör vi det i python,
    #   och kan de olika funktionerna ta emot olika många inputs?

    def snellslaw(self,r,s,n1,n2,n):
        if s.isinstance(Plane):
            pl = s
            # intersection between ray and plane
            p = pl.intersectionpoint(r)

            # normal line
            perp = pl.perpline(p)

            # Incidence angle
            theta_i = angle(r.toline(), perp)
            if theta_i > np.pi / 2:
                theta_i = np.pi - theta_i

            # Transmission angle
            theta_t = np.arcsin(n1 / n2 * np.sin(theta_i))

            # translation to origin(0)
            r0 = r.translate(-p)
            # pl0 = pl.translate(-p); %
            perp0 = perp.translate(-p)
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
            [Rs, Ts, rs, ts] = Ray.rtcoeffs(theta_i, n1, n2)

            # Rp and Tp(r3.pol.Vy and Vz)
            [Rp, Tp, rp, tp] = Ray.rtcoeffp(theta_i, n1, n2)

            # Reflected ray
            if r3.v.Z > 0:
                r_r3 = -r3.xrotation(2 * theta_i)
            else:
                r_r3 = -r3.xrotation(-2 * theta_i)

            r_r3.pol.Vx = abs(np.dot(rs,r_r3.pol.Vx))
            r_r3.pol.Vy = abs(np.dot(rp,r_r3.pol.Vy))
            r_r3.pol.Vz = abs(np.dot(rp,r_r3.pol.Vz))

            cs = np.divide((np.square(r3.pol.Vx)),np.square(r3.pol.norm()))
            cp = np.divide((np.square(r3.pol.Vy)) + np.square(r3.pol.Vz),np.square(r3.pol.norm()))
            R = np.dot(cs,Rs) + np.dot(cp,Rp)

            r_r3.v.X = np.zeros(np.shape(r))
            r_r3.v.Y = np.zeros(np.shape(r))
            r_r3.v.Z = np.zeros(np.shape(r))
            r_r3.P = np.dot(r.P,R)
            r_r3.pol.X = np.zeros(np.shape(r))
            r_r3.pol.Y = np.zeros(np.shape(r))
            r_r3.pol.Z = np.zeros(np.shape(r))

            # transmitted ray
            if r3.v.Z > 0:
                r_t3 = -r3.xrotation(-np.pi + theta_i - theta_t)
                r_t3.pol.Vy = -r_t3.pol.Vy
                r_t3.pol.Vz = -r_t3.pol.Vz
            else:
                r_t3 = -r3.xrotation(np.pi - theta_i + theta_t)
                r_t3.pol.Vy = -r_t3.pol.Vy
                r_t3.pol.Vz = -r_t3.pol.Vz

            r_t3.pol.Vx = np.dot(ts,r_t3.pol.Vx)
            r_t3.pol.Vy = np.dot(tp,r_t3.pol.Vy)
            r_t3.pol.Vz = np.dot(tp,r_t3.pol.Vz)

            T = 1 - R

            r_t3.v.X = np.zeros(np.shape(r))
            r_t3.v.Y = np.zeros(np.shape(r))
            r_t3.v.Z = np.zeros(np.shape(r))
            r_t3.P = np.dot(r.P,T)
            r_t3.pol.X = np.zeros(np.shape(r))
            r_t3.pol.Y = np.zeros(np.shape(r))
            r_t3.pol.Z = np.zeros(np.shape(r))

            # ManagesTIR
            indices = [1, 2, 3]
            tir = [j for (i, j) in zip(np.imag(theta_t), indices) if i != 0]
            tir = [x - 1 for x in tir]
            #tir = np.imag(theta_t) != 0

            r_r3.P(tir) = r.P(tir)
            r_r3.pol.Vx(tir) = np.abs(r_r3.pol.Vx(tir))
            r_r3.pol.Vy(tir) = np.abs(r_r3.pol.Vy(tir))
            r_r3.pol.Vz(tir) = np.abs(r_r3.pol.Vz(tir))

            r_t3.P(tir) = 0
            r_t3.v.Vx(tir) = np.real(r_t3.v.Vx(tir))
            r_t3.v.Vy(tir) = np.real(r_t3.v.Vy(tir))
            r_t3.v.Vz(tir) = np.real(r_t3.v.Vz(tir))
            r_t3.pol.Vx(tir) = np.abs(r_t3.pol.Vx(tir))
            r_t3.pol.Vy(tir) = np.abs(r_t3.pol.Vy(tir))
            r_t3.pol.Vz(tir) = np.abs(r_t3.pol.Vz(tir))

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

        elif s.isinstance(Superficies):
            # intersection between ray and sphere
            p = s.intersectionpoint(r, n)

            # normal line
            ln = s.perpline(p)

            # tangent plane
            pl = Plane.perpto(ln, p)

            #Snell's law
            output = r.snellslaw(pl, n1, n2)
            r_r = output['r_r']
            r_t = output['r_t']
            perp = output['perp']

        return {
            "r_r": r_r,
            "r_t": r_t,
            "perp": perp
        }

    @staticmethod
    def rtcoeffs(self,theta_i,n1,n2):
        theta_t = np.arcsin(n1 / n2 * np.sin(theta_i))

        rs = np.divide((n1 * np.cos(theta_i) - n2 * np.cos(theta_t)),(n1 * np.cos(theta_i) + n2 * np.cos(theta_t)))
        ts = np.divide(2 * n1 * np.cos(theta_i),(n1 * np.cos(theta_i) + n2 * np.cos(theta_t)))

        Rs = np.abs(np.square(rs))
        Ts = 1 - Rs

        return {
            "Rs": Rs,
            "Ts": Ts,
            "rs": rs,
            "ts": ts
        }

    @staticmethod
    def rtcoeffp(self,theta_i,n1,n2):
        theta_t = np.arcsin(n1 / n2 * np.sin(theta_i))

        rp = np.divide((n1 * np.cos(theta_t) - n2 * np.cos(theta_i)),(n1 * np.cos(theta_t) + n2 * np.cos(theta_i)))
        tp = np.divide(2 * n1 * np.cos(theta_i),(n1 * np.cos(theta_t) + n2 * np.cos(theta_i)))

        Rp = np.abs(np.square(rp))
        Tp = 1 - Rp

        return {
            "Rp": Rp,
            "Tp": Tp,
            "rp": rp,
            "tp": tp
        }

    @staticmethod
    def beam2rays(self,b):
#        if nargin < 3
#            z1 = 0;
#            z2 = 1e-6;
#        end

        X = np.dot(b.r,np.cos(b.phi))
        Y = b.r.dot(np.sin(b.phi))
        Z = 0 * b.r
        Vx = 0 * b.r
        Vy = 0 * b.r
        Vz = np.ones(np.shape(b.r))
        # TODO: Här vetefan vad jag ska skriva
        v = Vector([X;X], [Y;Y], [Z;Z], [Vx;Vx], [Vy;Vy], [Vz;Vz])

        [Ex, Ey] = Transform.Pol2CarVector(b.phi, b.Ephi, b.Er)
        # TODO: Här vetefan vad jag ska skriva
        pol = Vector([X;X], [Y;Y], [Z;Z], [real(Ex);imag(Ex)], [real(Ey);imag(Ey)], 0 * [b.r;b.r])

        dr = b.r(1, 2) - b.r(1, 1)
        dphi = b.phi(2, 1) - b.phi(1, 1)
        P = np.dot(b.intensity(),b.r) * dr * dphi
        # TODO: Här vetefan vad jag ska skriva
        P = [P. * np.abs(Ex). ^ 2. / (np.abs(Ex). ^ 2 + np.abs(Ey). ^ 2) P. * np.abs(Ey). ^ 2. / (np.abs(Ex). ^ 2 + np.abs(Ey). ^ 2)]

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

    def translate(self,r,dp):
        # TRANSLATE 3D translation of ray set
        #
        # Rt = TRANSLATE(R,dP) translates set of rays R by dP.
        #   If dP is a Point, the translation corresponds to the
        #   coordinates X, Y and Z.
        #   If dP is a Vector, the translation corresponds to the
        #   components Vx, Vy and Vz.
        #
        # See also Ray, Point, Vector.

        r_t = r
        r_t.v = r.v.translate(dp)
        r_t.pol = r.pol.translate(dp)
        return r_t

    def xrotation(self, r, phi):
        # XROTATION Rotation around x-axis of ray set
        #
        # Rr = XROTATION(R,phi) rotates set of rays R around x-axis
        #   by an angle phi [rad].
        #
        # See also Ray.

        r_r = r
        r_r.v = r.v.xrotation(phi)
        r_r.pol = r.pol.xrotation(phi)
        return r_r

    def yrotation(self, r, phi):
        # YROTATION Rotation around y-axis of ray set
        #
        # Rr = YROTATION(R,phi) rotates set of rays R around y-axis
        #   by an angle phi [rad].
        #
        # See also Ray.

        r_r = r
        r_r.v = r.v.yrotation(phi)
        r_r.pol = r.pol.yrotation(phi)
        return r_r

    def zrotation(self, r, phi):
        # ZROTATION Rotation around z-axis of ray set
        #
        # Rr = ZROTATION(R,phi) rotates set of rays R around z-axis
        #   by an angle phi [rad].
        #
        # See also Ray.

        r_r = r
        r_r.v = r.v.zrotation(phi)
        r_r.pol = r.pol.zrotation(phi)
        return r_r

    def numel(self, r):
        # NUMEL Number of rays
        #
        # N = NUMEL(R) number of rays in set R.
        #
        # See also Ray.

        return r.v.size

    def size(self, r, varargin):
        # SIZE Size of the ray set
        #
        # S = SIZE(R) returns a two-element row vector with the number
        #   of rows and columns in the ray set R.
        #
        # S = SIZE(R,DIM) returns the length of the dimension specified
        #   by the scalar DIM in the ray set R.
        #
        # See also Ray.

        if varargin.size > 0:
            s = r.v.size(varargin[1])
        else:
            s = r.v.size()
        return s

    def uminus(self, r):
        # UMINUS Unitary minus (components)
        #
        # Rm = UMINUS(R) Unitary minus (Rm = -R).
        #   -R inverts the components of R.
        #   The points of application and polarizations are left unchanged.
        #
        # See also Ray.

        r_m = r
        r_m.v = -r_m.v
        return r_m

    def angle(self, r1, r2):
        # ANGLE Angle (coordinates)
        #
        # PHI = ANGLE(R1,R2) calculates the angle between the set of
        #   rays R1 and R2.
        #
        # See also Ray.

        return angle(r1.v, r2.v)

    def versor(self, r):
        # VERSOR Unitary vector
        #
        # U = VERSOR(R) returns the unit vector set corresponding to
        #   the ray set R.
        #   The coordinates of U are the points of application of R.
        #
        # See also Ray, Vector.

        return r.v.versor()

    def toline(self, r):
        # TOLINE Ray to line
        #
        # LN = TOLINE(R) converts the set of rays R into the set of
        #   lines LN. The coordinates X, Y and Z of the initial points
        #   of the lines are the points of application of R and
        #   the coordinates of the final points are the sum of the coordiantes
        #   of the initial points and of the components of R.
        #
        # See also Ray, SLine.

        return r.v.toline()