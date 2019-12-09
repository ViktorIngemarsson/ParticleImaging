import math

def dotproduct(v1, v2):
    return sum((a * b) for a, b in zip(v1, v2))

def length(v):
    return math.sqrt(dotproduct(v, v))


def angle(v1, v2):
    return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))


class Ray:
    def __init__(self, vectors, power):
        self.vectors = vectors
        self.power = power

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

    def snellslaw(self, r,s,n1,n2,n):
        # SNELLSLAW Snell's law: reflected and transmitted rays at a surface
        #
        # [Rr,Rt,PERP] = SNELLSLAW(R,S,n1,n2) calculates the reflected
        #   ray Rr, and the transmitted ray set Rt and the
        #   perpendicular line set PERP for a ray set R incident to a
        #   superficies S. n1 and n2 represents the refractive index
        #   in the incoming medium and in the transmission medium.
        #
        # [Rr,Rt,PERP] = SNELLSLAW(R,S,n1,n2,n) n [default = 1]
        #   defines what intersaction between the ray and the surface
        #   should be used.
        #
        # See also Ray, Superficies, Plane, Spherical, Ellipsoidal, Cylindrical.

        if nargin < 5
            n = 1;
        end