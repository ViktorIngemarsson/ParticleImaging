from Point import Point
from Vector import Vector
from Ray import Ray
from ParticleSpherical import ParticleSpherical
import numpy as np

# Parameters
# Medium
nm = 1.33  # Medium refractive index

# Particle
R = 1e-6  # Particle radius [m]
nP = 1.50  # Particle refractive index
c = Point(0, 0, 0)  # Particle center [m]

# Initialization
# Particle
bead = ParticleSpherical(c, R, nm, nP)

theta = 0.5 #% [0: 1:89.9]./ 180 * np.pi

# Ray
v = Vector(-2 * R, R * np.sin(theta), 0, 1, 0.2, 0.2) # Direction
P = 1 # Power[W]
pol = Vector(0, 0, 0, 0, 1, 1)
pol = v.mtimes(pol)
pol = pol.versor() # Polarization
r = Ray(v, P, pol)

# Scattered rays
s = bead.scattering(r, 1e-12, 6)

# Scattering coefficients
#f = bead.force(r, 1e-12, 6);

# #Plotting
vectorsTransmitted = list()
vectorsReflected = list()
for vector in s:
    vectorsTransmitted.append(vector["r_t"].v)
    vectorsReflected.append(vector["r_r"].v)

s_t = np.asarray(vectorsTransmitted)
s_r = np.asarray(vectorsReflected)
s_t[0].plot_multiple_vectors(s_t,R)
s_r[0].plot_multiple_vectors(s_r,R)

