import numpy as np

def TranslateCoordinates(x,y,z,xRotation,zRotation):
    '''
    Rotates coordinates round axis, x and y, with defines degrees...
    xRotation and zRotation.
    :param x: coodinate x
    :param y: coodinate y
    :param z: coodinate z
    :param xRotation: rotation in dg around x-axis
    :param zRotation: rotation in dg around z-axis
    :return: coordinates in new coordinate system
    '''
    xPrim = np.cos(zRotation)*x + np.sin(zRotation)*y
    yPrim = -np.sin(zRotation)*x + np.cos(zRotation)*y
    yPrim = np.cos(xRotation)*yPrim - np.sin(xRotation)*z
    zPrim = np.sin(xRotation)*yPrim + np.cos(xRotation)*z
    return xPrim, yPrim, zPrim