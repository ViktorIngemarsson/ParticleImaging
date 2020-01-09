from time import process_time
from GeneratingOneImage import GeneratingOneImage
from TranslateCoordinates import TranslateCoordinates
from generatingPerlinNoise import finalNoiseGeneration
import matplotlib.pyplot as plt

################################################## Parameters ##########################################################

# Medium
nm = 1.33  # Medium refractive index

# Particle
R = 0.000005  # Particle radius [m]
nP = 1.50  # Particle refractive index
rho = 360 * 10e10

# Particle center [m]
# The particle center is for now not changeable
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

# Light
scattering_max_number_of_iterations = 6

# Camera
lens_size_x = 0.001  # 1mm
lens_size_y = 0.001
resolution = 1000  # Assuming quadratic image
n_of_cameras = 3
xRotation = [0, 120, 30]
yRotation = [0, 0, 115]

################################################## Generate images #####################################################
# Generate images and plot them.
image = []
f, axarr = plt.subplots(1, 3)
previous_coordinates = [particle_center_x, particle_center_y, particle_center_z]

# Images from the cameras
for camera in range(n_of_cameras):
    t1_start = process_time() # Timing the program for one image
    [new_x, new_y, new_z] = TranslateCoordinates(previous_coordinates[0], previous_coordinates[1],
                                                 previous_coordinates[2], xRotation[camera], yRotation[camera])
    image.append(GeneratingOneImage(new_x, new_y, R, rho, nm, nP, pol_X, pol_Y, pol_Z, pol_Vx, pol_Vy, pol_Vz,
                                      resolution, lens_size_x, lens_size_y, scattering_max_number_of_iterations) + \
                   finalNoiseGeneration(resolution))
    # Timer
    t1_stop = process_time()
    print("Elapsed CPU time (s) to generate image {}:".format(camera+1), t1_stop - t1_start)

    previous_coordinates = [new_x, new_y, new_z]
    axarr[camera].imshow(image[camera], cmap='gray_r', vmin=0, vmax=1)
    axarr[camera].axis('off')
    axarr[camera].title.set_text('Camera{} ' .format(camera+1))

plt.show()
