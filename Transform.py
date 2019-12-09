import numpy as np

# TODO: Add more transforms.

class Transform:
    def Pol2CarVector(self, phi,Vphi,Vr):
        return {
        "Vx": np.dot(np.cos(phi),Vr) - np.dot(np.sin(phi),Vphi),
        "Vy": np.dot(np.sin(phi),Vr) + np.dot(np.cos(phi),Vphi)
        }