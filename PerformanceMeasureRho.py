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
rho = [1e11, 1e12, 1e13, 1e14, 1e15, 1e16]

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
resolution = 600

# Number of cameras to capture images of droplet.
n_of_cameras = 3

# The rotation required to flip the coordinates from one camera to the other.
xRotation = [0, 120, 30]
yRotation = [0, 0, 115]

# Generate images and plot them.
image = []
#f, axarr = plt.subplots(1, 3)
previous_coordinates = [particle_center_x, particle_center_y, particle_center_z]

reruns = 5
timeaverage = 0
timeAverageTotal = np.zeros(10)
indexing = 0
# Images from the cameras
for iRho in rho:
    for iRerun in range(reruns):
        for camera in range(n_of_cameras):
            t1_start = process_time()  # Timing the program for one image
            [new_x, new_y, new_z] = TranslateCoordinates(previous_coordinates[0], previous_coordinates[1],
                                                 previous_coordinates[2], xRotation[camera], yRotation[camera])
            image.append(GeneratingOneImage(new_x, new_y, R, iRho, nm, nP, pol_X, pol_Y, pol_Z, pol_Vx, pol_Vy, pol_Vz,
                                    resolution, lens_size, scattering_max_number_of_iterations) + \
                 FinalNoiseGeneration(resolution, noise_intensity))
    # Timer
            t1_stop = process_time()
            timeaverage = timeaverage +  t1_stop - t1_start
    #print("Elapsed CPU time (s) to generate image {}:".format(camera + 1), t1_stop - t1_start)

            previous_coordinates = [new_x, new_y, new_z]
   # axarr[camera].imshow(image[camera], cmap='gray_r', vmin=0, vmax=1)
  #  axarr[camera].axis('off')
 #   axarr[camera].title.set_text('Camera{} '.format(camera + 1))
    timeAverageTotal[indexing] = timeaverage/reruns
    indexing = indexing+1
    timeaverage = 0
#plt.show()
print(timeAverageTotal)
plt.plot(rho, timeAverageTotal[0:6])

# naming the x axis
plt.xlabel('Rays per sqm')
# naming the y axis
plt.ylabel('Time to generate one observation [s]')

# function to show the plot
plt.show()