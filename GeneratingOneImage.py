from Point import Point
from Vector import Vector
from Ray import Ray
from ParticleSpherical import ParticleSpherical
import numpy as np
from LinePlaneCollision import LinePlaneCollision
from GeneratingPoints import GeneratingCoordinates


def GeneratingOneImage(particle_center_x, particle_center_y, R, rho, nm, nP, pol_X, pol_Y, pol_Z,
                       pol_Vx, pol_Vy, pol_Vz, resolution, lens_size, scattering_number_of_iterations):

    """
    This method generates an image based on the rays propagated from the light source, refracted through a particle,
        and towards a camera lens.
    :param particle_center_x: Particle center x coordinate.
    :param particle_center_y: Particle center y coordinate.
    :param R: Particle radius [m].
    :param rho: Ray density per square meter.
    :param nm: Refractive index of the medium before ray impact.
    :param nP: Refractive index of the medium after ray impact.
    :param pol_X: Polarization
    :param pol_Y: Polarization
    :param pol_Z: Polarization
    :param pol_Vx: Polarization
    :param pol_Vy: Polarization
    :param pol_Vz: Polarization
    :param resolution: Resolution of the generated image. Assuming quadratic image.
    :param lens_size: Camera lens size [m], assuming quadratic lens.
    :param scattering_number_of_iterations: # Maximum number of iterations the rays are scattered within the particle.
    :return: The particle captures by the camera as an image.
    """

    # Produce one image.
    c = Point(0, 0, 0)
    bead = ParticleSpherical(c, R, nm, nP)

    pixels = np.zeros(resolution * resolution)

    # Coordinates for the rays hitting the sphere
    x, y = GeneratingCoordinates(rho, R, 0, 0)
    number_of_rays = np.size(x)
    P = np.ones(number_of_rays)  # Power[W]
    z = np.zeros(number_of_rays)
    Vx = np.zeros(number_of_rays)
    Vy = np.zeros(number_of_rays)
    Vz = np.ones(number_of_rays)

    v = Vector(x * R, y * R, z, Vx, Vy, Vz)

    pol = Vector(pol_X, pol_Y, pol_Z, pol_Vx, pol_Vy, pol_Vz)
    pol = v.mtimes(pol)
    pol = pol.versor()

    r = Ray(v, P, pol)
    s = bead.scattering(r=r, N=scattering_number_of_iterations)

    plane_normal = np.array([0, 0, 1])  # Define lens plane
    plane_point = np.array([0, 0, 1])  # Any point on the plane

    exiting_rays = s[-1]["t"]
    ray_power = exiting_rays.P

    for i in range(0, np.size(exiting_rays.v.Vx) - 1):
        ray_direction = np.array([exiting_rays.v.Vx[i], exiting_rays.v.Vy[i], exiting_rays.v.Vz[i]])
        ray_point = np.array([exiting_rays.v.X[i], exiting_rays.v.Y[i], exiting_rays.v.Z[i]]) + np.array(
            [particle_center_x, particle_center_y, 0])

        psi = LinePlaneCollision(plane_normal, plane_point, ray_direction, ray_point, lens_size)
        if psi is not None:  # Ray hits the camera lens
            x_pixel = np.round((psi[0] + lens_size / 2) / lens_size * (resolution-1))
            y_pixel = np.round((psi[1] + lens_size / 2) / lens_size * resolution)
            pixel_index = y_pixel * resolution + x_pixel
            pixels[int(pixel_index)] = pixels[int(pixel_index)] + ray_power[i]

    # Normalize the sphere-rays
    if max(pixels) > 0:
        pixels = pixels / max(pixels)

    return np.reshape(pixels, (resolution, resolution))
