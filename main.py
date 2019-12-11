from Point import Point
from Vector import Vector
from Ray import Ray
from ParticleSpherical import ParticleSpherical
import numpy as np
from LinePlaneCollision import LinePlaneCollision

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

#theta = 0.5 #% [0: 1:89.9]./ 180 * np.pi

# Ray
#v = Vector(-2 * R, R * np.sin(theta), 0, 1, 0.2, 0.2) # Direction
#v = Vector(np.asarray([R/2,R/3]),np.asarray([R/2,R/3]),np.asarray([0,0]),np.asarray([0,0]),np.asarray([0,0]),np.asarray([1,1]))
v = Vector(0,0,0,0,0,1)
P = 1 # Power[W]
pol = Vector(0, 0, 0, 0, 1, 1)
pol = v.mtimes(pol)
pol = pol.versor() # Polarization
r = Ray(v, P, pol)

# Scattered rays
s = bead.scattering(r, 1e-12, 6)

# Scattering coefficients
#f = bead.force(r, 1e-12, 6);

## Plotting one ray with a sphere
vectorsTransmitted = list()
vectorsReflected = list()
for vector in s:
    vectorsTransmitted.append(vector["r_t"].v)
    vectorsReflected.append(vector["r_r"].v)

s_t = np.asarray(vectorsTransmitted)
s_r = np.asarray(vectorsReflected)
s_t[0].plot_multiple_vectors(s_t,R)
s_r[0].plot_multiple_vectors(s_r,R)

## Does the rays hit the camera-lens?
resolution = 720
lens_size_x = 0.01
lens_size_y = 0.01

# Define plane
planeNormal = np.array([0, 0, 1])
planePoint = np.array([0, 0, 1])  # Any point on the plane

for vector in s:
    rayDirection = np.array([vector["r_r"].v.Vx, vector["r_r"].v.Vy, vector["r_r"].v.Vz])
    rayPoint = np.array([vector["r_r"].v.X, vector["r_r"].v.Y, vector["r_r"].v.Z])

    Psi = LinePlaneCollision([planeNormal], planePoint, rayDirection, rayPoint, lens_size_x, lens_size_y)
    if Psi is None:
        print("Intersection at", Psi,", Missed the camera!")
    else:
        print("Intersection at", Psi, "Hit the camera!")
        pixels = np.zeros(resolution * resolution)
        x_interval_length = lens_size_x / resolution
        y_interval_length = lens_size_y / resolution
        x_pixel = np.round((Psi[0] + lens_size_x / 2) / lens_size_x * resolution)
        y_pixel = np.round((Psi[1] + lens_size_y / 2) / lens_size_y * resolution)
        pixel_index = y_pixel * resolution + x_pixel
        pixels[int(pixel_index)] = + 1






