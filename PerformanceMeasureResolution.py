from time import process_time
from GeneratingOneImage import GeneratingOneImage
from TranslateCoordinates import TranslateCoordinates
from GeneratingNoise import FinalNoiseGeneration
import matplotlib.pyplot as plt
import numpy as np

"""
Example script that generates three images, from three cameras, and plots them. Parameters are explained and specified
 below. 
"""

# Refractive index of the medium around the particle (ie. air).
nm = 1.33

# Particle radius [m].
R = 0.000002

# Refractive index of particle medium (ie. water).
nP = 1.50

# Ray density per square meter.
rho = 1e14

# Particle center coordinates.
particle_center_x = 0.00014
particle_center_y = 0.00011
particle_center_z = 0

# Polarization
pol_X = 0
pol_Y = 0
pol_Z = 0
pol_Vx = 0
pol_Vy = 1
pol_Vz = 1

# Maximum number of iterations the rays are scattered within the particle.
scattering_max_number_of_iterations = 6

# Perlin noise with an intensity is generated.
noise_intensity = 0.2

# Camera lens size [m]. This is equal to the image width. Also assuming quadratic lens/image.
lens_size = 0.001

# Resolution of the generated image. Assuming quadratic image. Number of pixels in image becomes resolution*resolution.
resolution = [100, 200, 400, 800, 1600, 3200]
reruns = 5

# Number of cameras to capture images of droplet.
n_of_cameras = 3

# The rotation required to flip the coordinates from one camera to the other.
xRotation = [0, 120, 30]
yRotation = [0, 0, 115]

# Generate images and plot them.
image = []
#f, axarr = plt.subplots(1, 3)
previous_coordinates = [particle_center_x, particle_center_y, particle_center_z]
indexing = 0
totalTimeAverage = np.zeros(7)
timeAverage = 0
# Images from the cameras
for iResolution in resolution:
    indexing = indexing + 1
    for iRerun in range(reruns):
        for camera in range(n_of_cameras):
            t1_start = process_time()  # Timing the program for one image
            [new_x, new_y, new_z] = TranslateCoordinates(previous_coordinates[0], previous_coordinates[1],
                                                 previous_coordinates[2], xRotation[camera], yRotation[camera])
            image.append(GeneratingOneImage(new_x, new_y, R, rho, nm, nP, pol_X, pol_Y, pol_Z, pol_Vx, pol_Vy, pol_Vz,
                                    iResolution, lens_size, scattering_max_number_of_iterations) + \
                 FinalNoiseGeneration(iResolution, noise_intensity))
            t1_stop = process_time()
            timeAverage = timeAverage + (t1_stop - t1_start)/3

            previous_coordinates = [new_x, new_y, new_z]

    totalTimeAverage[indexing] = timeAverage/reruns
    timeAverage = 0

print(totalTimeAverage)

plt.plot(resolution, totalTimeAverage[1:7])
plt.xlabel('Resolution')
plt.ylabel('Time to generate one observation [s]')
plt.show()

