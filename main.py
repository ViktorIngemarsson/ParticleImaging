from Point import Point
from Vector import Vector
from Ray import Ray
from ParticleSpherical import ParticleSpherical
import numpy as np
from LinePlaneCollision import LinePlaneCollision
from GeneratingPoints import GeneratingCoordinates
import matplotlib.pyplot as plt
from time import process_time
from GeneratingOneImage import GeneratingOneImage
from TranslateCoordinates import TranslateCoordinates
from generatingPerlinNoise import finalNoiseGeneration
import matplotlib.cm as cm

################################################## Parameters ##########################################################

# Medium
nm = 1.33  # Medium refractive index

# Particle
R = 0.00005  # Particle radius [m]
nP = 1.50  # Particle refractive index
rho = 360 * 10e9

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
scattering_number_of_iterations = 6

# Camera
lens_size_x = 0.001  # 1mm
lens_size_y = 0.001
resolution = 400  # Assuming quadratic image

# Camera 2 relative camera 1
xRotation = 90
yRotation = 0


###################################################
#Timing the program
###################################################
t1_start = process_time()
###################################################
GeneratingOneImage(particle_center_x, particle_center_y, particle_center_z, R,rho, nm, nP,pol_X, pol_Y, pol_Z, pol_Vx, pol_Vy, pol_Vz,resolution,lens_size_x,lens_size_y,scattering_number_of_iterations)
###################################################
#Timing the program
###################################################
t1_stop = process_time()
print("Elapsed CPU time to generate one image:", t1_stop-t1_start)
###################################################

[a,b,c] = TranslateCoordinates(particle_center_x,particle_center_y,particle_center_z,xRotation,yRotation)
print(a)
print(b)
print(c)
img = finalNoiseGeneration(1000)
plt.imshow(img, cmap=cm.gray)
plt.show()
########################################################################################################################

