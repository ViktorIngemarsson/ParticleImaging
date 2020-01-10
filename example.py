from time import process_time
from GeneratingOneImage import GeneratingOneImage
from TranslateCoordinates import TranslateCoordinates
from GeneratingPerlinNoise import FinalNoiseGeneration
import matplotlib.pyplot as plt

"""
Example script that generates three images, from three cameras, and plots them. Parameters are explained and specified
 below. 
"""

# Refractive index of the medium around the particle (ie. air).
n_m = 1.0

# Refractive index of particle medium (ie. water).
n_p = 1.35

# Particle radius [m].
r = 0.000002

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
resolution = 400

# Number of cameras to capture images of droplet.
n_of_cameras = 3

# The rotation required to flip the coordinates from one camera to the other.
x_rotation = [0, 120, 30]
y_rotation = [0, 0, 115]

# Generate images and plot them.
image = []
f, axarr = plt.subplots(1, 3)
previous_coordinates = [particle_center_x, particle_center_y, particle_center_z]

# Images from the cameras
for camera in range(n_of_cameras):
    t1_start = process_time()  # Timing the program for one image
    [new_x, new_y, new_z] = TranslateCoordinates(previous_coordinates[0], previous_coordinates[1],
                                                 previous_coordinates[2], x_rotation[camera], y_rotation[camera])
    image.append(GeneratingOneImage(new_x, new_y, r, rho, n_m, n_p, pol_X, pol_Y, pol_Z, pol_Vx, pol_Vy, pol_Vz,
                                    resolution, lens_size, scattering_max_number_of_iterations) + \
                 FinalNoiseGeneration(resolution, noise_intensity))
    # Timer
    t1_stop = process_time()
    print("Elapsed CPU time (s) to generate image {}:".format(camera + 1), t1_stop - t1_start)

    previous_coordinates = [new_x, new_y, new_z]
    axarr[camera].imshow(image[camera], cmap='gray_r', vmin=0, vmax=1)
    axarr[camera].axis('off')
    axarr[camera].title.set_text('Camera{} '.format(camera + 1))

plt.show()

