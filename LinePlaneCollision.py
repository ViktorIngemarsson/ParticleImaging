from __future__ import print_function
import numpy as np


def LinePlaneCollision(planeNormal, planePoint, rayDirection, rayPoint, planeSizeX, planeSizeY, epsilon=1e-6):
    ndotu = np.dot(planeNormal, rayDirection)
    if np.abs(ndotu) < epsilon:
        Psi = None
    else:
        w = rayPoint - planePoint
        si = -np.dot(planeNormal, w) / ndotu
        Psi = w + si * rayDirection + planePoint

    #Allowed limits=
    xMin = -planeSizeX/2
    xMax = planeSizeX/2
    yMin = -planeSizeY/2
    yMax = planeSizeY/2

    if Psi[0] < xMin or Psi[0] > xMax or Psi[1] < yMin or Psi[1] > yMax:
        Psi = None

    return Psi