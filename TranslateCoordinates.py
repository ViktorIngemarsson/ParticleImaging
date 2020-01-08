import numpy as np

def TranslateCoordinates(x,y,z,xRotation,zRotation):
    xPrim = np.cos(zRotation)*x + np.sin(zRotation)*y
    yPrim = -np.sin(zRotation)*x + np.cos(zRotation)*y
    yPrim = np.cos(xRotation)*yPrim - np.sin(xRotation)*z
    zPrim = np.sin(xRotation)*yPrim + np.cos(xRotation)*z
    return xPrim, yPrim, zPrim