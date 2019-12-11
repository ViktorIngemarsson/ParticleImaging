from Point import Point
from Vector import Vector
from Ray import Ray
from ParticleSpherical import ParticleSpherical
import numpy as np
from LinePlaneCollision import LinePlaneCollision
from GeneratingPoints import GeneratingCoordinates
import matplotlib.pyplot as plt

# Parameters
# Medium
nm = 1.33  # Medium refractive index

# Particle
R = 0.01  # Particle radius [m]
nP = 1.50  # Particle refractive index
c = Point(0, 0, 0)  # Particle center [m]

# Initialization
# Particle
bead = ParticleSpherical(c, R, nm, nP)
P = 1  # Power[W]
pol = Vector(0, 0, 0, 0, 1, 1)

number_of_rays = 500
resolution = 100
pixels = np.zeros(resolution * resolution)
x, y = GeneratingCoordinates(number_of_rays)
x = x * R
y = y * R

lens_size_x = 0.001
lens_size_y = 0.001

# Define lens plane
planeNormal = np.array([0, 0, 1])
planePoint = np.array([0, 0, 1])  # Any point on the plane

for i in range(0, number_of_rays - 1):
    print(i)
    v = Vector(x[i], y[i], 0, 0, 0, 1)

    pol = v.mtimes(pol)
    pol = pol.versor()  # Polarization
    r = Ray(v, P, pol)

    # Scattered rays
    s = bead.scattering(r, 1e-12, 6)

    for vector in s:
        rayDirection = np.array([vector["r_r"].v.Vx, vector["r_r"].v.Vy, vector["r_r"].v.Vz])
        rayPoint = np.array([vector["r_r"].v.X, vector["r_r"].v.Y, vector["r_r"].v.Z])

        Psi = LinePlaneCollision([planeNormal], planePoint, rayDirection, rayPoint, lens_size_x, lens_size_y)
        if Psi is not None: #Ray hits the camera lens
            x_interval_length = lens_size_x / resolution
            y_interval_length = lens_size_y / resolution
            x_pixel = np.round((Psi[0] + lens_size_x / 2) / lens_size_x * resolution)
            y_pixel = np.round((Psi[1] + lens_size_y / 2) / lens_size_y * resolution)
            pixel_index = y_pixel * resolution + x_pixel
            pixels[int(pixel_index)] = pixels[int(pixel_index)] + 1

norm_int = max(pixels)
pixels = pixels / norm_int
image = np.reshape(pixels, (resolution, resolution))

plt.imshow(image, cmap='gray_r', vmin=0, vmax=1)
plt.show()
