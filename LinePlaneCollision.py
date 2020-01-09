from __future__ import print_function
import numpy as np


def LinePlaneCollision(plane_normal, plane_point, ray_direction, ray_point, plane_size, epsilon=1e-6):

    """
    :param plane_normal: Normal vector to the plane
    :param plane_point: A point in the plane
    :param ray_direction: The direction of the ray
    :param ray_point: A point along the ray
    :param plane_size: The size of the quadratic plane
    :param epsilon:
    :return: Coordinates of the interception point. Returns None if the ray misses the plane_size.
    """

    ndotu = np.dot(plane_normal, ray_direction)
    if np.abs(ndotu) < epsilon:
        psi = None
    else:
        w = ray_point - plane_point
        si = -np.dot(plane_normal, w) / ndotu
        psi = w + si * ray_direction + plane_point

    x_min = -plane_size/2
    x_max = plane_size/2
    y_min = -plane_size/2
    y_max = plane_size/2

    if psi[0] < x_min or psi[0] > x_max or psi[1] < y_min or psi[1] > y_max:
        psi = None

    return psi