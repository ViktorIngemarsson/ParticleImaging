from Point import Point
from Vector import Vector
from Ray import Ray
from ParticleSpherical import ParticleSpherical
import numpy as np
from LinePlaneCollision import LinePlaneCollision
from GeneratingPoints import GeneratingCoordinates
import matplotlib.pyplot as plt


def GeneratingOneImage(particle_center_x, particle_center_y, particle_center_z, R,rho, nm, nP,pol_X, pol_Y, pol_Z, pol_Vx, pol_Vy, pol_Vz,resolution,lens_size_x,lens_size_y,scattering_number_of_iterations):
    ### Produce one image.
    c = Point(0, 0, 0)
    bead = ParticleSpherical(c, R, nm, nP)

    pixels = np.zeros(resolution * resolution)

    # Coordinates for the rays hitting the sphere
    testX, testY = GeneratingCoordinates(rho, R, 0, 0)
    number_of_rays = np.size(testX)
    P = np.ones(number_of_rays)  # Power[W]
    testZ = np.zeros(number_of_rays)
    testVx = np.zeros(number_of_rays)
    testVy = np.zeros(number_of_rays)
    testVz = np.ones(number_of_rays)

    # TODO: Här borde man inte multiplicera med R ....
    v = Vector(testX * R, testY * R, testZ, testVx, testVy, testVz)  # light source is 0.5m from sphere (z)

    pol = Vector(pol_X, pol_Y, pol_Z, pol_Vx, pol_Vy, pol_Vz)
    pol = v.mtimes(pol)
    pol = pol.versor()  # Polarization

    r = Ray(v, P, pol)
    s = bead.scattering(r, scattering_number_of_iterations)

    # Define lens plane
    planeNormal = np.array([0, 0, 1])
    planePoint = np.array([0, 0, 1])  # Any point on the plane

    exiting_rays = s[-1]["t"]
    rayPowers = exiting_rays.P

    for i in range(0, np.size(exiting_rays.v.Vx) - 1):
        print("Computing trajectory for ray:", i)
        rayDirection = np.array([exiting_rays.v.Vx[i], exiting_rays.v.Vy[i], exiting_rays.v.Vz[i]])
        rayPoint = np.array([exiting_rays.v.X[i], exiting_rays.v.Y[i], exiting_rays.v.Z[i]]) + np.array([particle_center_x,particle_center_y,0])

        Psi = LinePlaneCollision(planeNormal, planePoint, rayDirection, rayPoint, lens_size_x, lens_size_y)
        if Psi is not None:  # Ray hits the camera lens
            x_pixel = np.round((Psi[0] + lens_size_x / 2) / lens_size_x * (resolution-1))
            y_pixel = np.round((Psi[1] + lens_size_y / 2) / lens_size_y * (resolution))
            pixel_index = (y_pixel) * resolution + x_pixel
            pixels[int(pixel_index)] = pixels[int(pixel_index)] + rayPowers[i]

    # Normalize the sphere-rays
    if max(pixels) > 0:
        pixels = pixels / max(pixels)

    image = np.reshape(pixels, (resolution, resolution))
    fig = plt.figure(figsize=(1, 3))
    plt.imshow(image, cmap='gray_r', vmin=0, vmax=1)
    fig.suptitle('Ray optics for the sphere', fontsize=20)
    plt.axis('off')
    plt.show()