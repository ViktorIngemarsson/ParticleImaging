from Point import Point
from Vector import Vector
from Ray import Ray
from ParticleSpherical import ParticleSpherical
import numpy as np
from LinePlaneCollision import LinePlaneCollision
from GeneratingPoints import GeneratingCoordinates
import matplotlib.pyplot as plt
from time import process_time

################################################## Parameters ##########################################################

# Medium
nm = 1.33  # Medium refractive index

# Particle
R = 0.00001  # Particle radius [m]
nP = 1.50  # Particle refractive index

# Particle center [m]
# The particle center is for now not changeable
particle_center_x = 0
particle_center_y = 0
particle_center_z = 0

# Polarization
pol_X = 0
pol_Y = 0
pol_Z = 0
pol_Vx = 0
pol_Vy = 1
pol_Vz = 1

# Light
number_of_rays = 100  # Number of rays that hits the sphere.
P = 1  # Power[W]
scattering_number_of_iterations = 6

# Camera
lens_size_x = 0.001  # 1mm
lens_size_y = 0.001
resolution = 1000  # Assuming quadratic image

########################################################################################################################

### Produce one image.


###################################################
#Timing the program
###################################################
t1_start = process_time()
###################################################

# Initialization Particle
c = Point(particle_center_x, particle_center_y, particle_center_z)
bead = ParticleSpherical(c, R, nm, nP)

# Polarization
pol = Vector(pol_X, pol_Y, pol_Z, pol_Vx, pol_Vy, pol_Vz)

pixels = np.zeros(resolution * resolution)
pixels_background = np.zeros(resolution * resolution)

# Coordinates for the rays hitting the sphere
rho = 360*10e9
x, y = GeneratingCoordinates(rho, R, particle_center_x, particle_center_y)
number_of_rays = np.size(x)
# Define lens plane
planeNormal = np.array([0, 0, 1])
planePoint = np.array([0, 0, 1])  # Any point on the plane
x_interval_length = lens_size_x / resolution
y_interval_length = lens_size_y / resolution

for i in range(0, number_of_rays - 1):
    print("Computing trajectory for ray:", i)

    # It is unknown why x and y have to be multiplied by R here, will be fixed in the future.
    v = Vector(x[i] * R, y[i] * R, -0.5, 0, 0, 1)  # light source is 0.5m from sphere (z)

    pol = v.mtimes(pol)
    pol = pol.versor()  # Polarization
    r = Ray(v, P, pol)

    # Scattered rays
    s = bead.scattering(r, scattering_number_of_iterations)

    for vector in s:
        rayDirection = np.array([vector["r_r"].v.Vx, vector["r_r"].v.Vy, vector["r_r"].v.Vz])
        rayPoint = np.array([vector["r_r"].v.X, vector["r_r"].v.Y, vector["r_r"].v.Z])

        Psi = LinePlaneCollision([planeNormal], planePoint, rayDirection, rayPoint, lens_size_x, lens_size_y)
        if Psi is not None:  # Ray hits the camera lens
            x_pixel = np.round((Psi[0] + lens_size_x / 2) / lens_size_x * resolution)
            y_pixel = np.round((Psi[1] + lens_size_y / 2) / lens_size_y * resolution)
            pixel_index = (y_pixel - 1) * resolution + x_pixel
            pixels[int(pixel_index)] = pixels[int(pixel_index)] + vector["r_r"].P

    x_pixel_background = np.round((x[i] + lens_size_x / 2) / lens_size_x * resolution)
    y_pixel_background = np.round((y[i] + lens_size_y / 2) / lens_size_y * resolution)
    pixel_index_background = (y_pixel_background - 1) * resolution + x_pixel_background
    pixels_background[int(pixel_index_background)] = pixels_background[int(pixel_index_background)] + 1

# Make light background coherent
for i in range(0, pixels_background.size - 1):
    if pixels_background[i] > 0:
        pixels_background[i] = 1
pixels_background = 1 - pixels_background

# Normalize the sphere-rays
pixels = pixels / max(pixels)

# Add together the sphere-rays with all other light
new_pixels = np.add(pixels_background, pixels)

# Normalize the image
new_pixels = new_pixels / max(new_pixels)

image = np.reshape(pixels, (resolution, resolution))

###################################################
#Timing the program
###################################################
t1_stop = process_time()
print("Elapsed CPU time to generate one image:", t1_stop-t1_start)
###################################################


fig = plt.figure(figsize=(1, 3))
plt.imshow(image, cmap='gray_r', vmin=0, vmax=1)
fig.suptitle('Ray optics for the sphere', fontsize=20)
plt.axis('off')
plt.show()
