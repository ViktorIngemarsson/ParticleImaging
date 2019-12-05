import numpy as np

class Droplet:
    def __init__(self, r, pos):
        self.radius = r
        self.pos = pos

class Ray:
    def __init__(self, pos, vectors, powers):
        self.outGoingPos = pos
        self.vector = vectors
        self.power = powers

    def inflectionWithDroplet(self, droplet):
        a = np.multiply(self.vector.v,self.vector.v)
        b = np.dot(np.dot(2,(np.subtract(self.outGoingPos.v, droplet.pos))),self.vector.v)
        c = np.subtract(np.dot(np.subtract(self.outGoingPos.v,droplet.pos),np.subtract(self.outGoingPos.v,droplet.pos)),droplet.radius*droplet.radius)
        delta = np.subtract(np.dot(b,b),np.dot(4,np.dot(a,c)))
        solution1 =  np.subtract(-b,np.divide(np.sqrt(delta),np.dot(2,a)))
        solution2 =  np.add(-b,np.divide(np.sqrt(delta),np.dot(2,a)))
        output1 = np.add(self.outGoingPos.v,np.multiply(solution1,self.vector.v))
        output2 = np.add(self.outGoingPos.v,np.multiply(solution2,self.vector.v))

    def getSphereNorm(self, droplet):
        sphereNorm = np.array([np.divide(np.subtract(self.outGoingPos[:,0], droplet.pos[0]), droplet.radius), np.divide(np.subtract(self.outGoingPos[:,1], droplet.pos[1]), droplet.radius), np.divide(np.subtract(self.outGoingPos[:,2], droplet.pos[2]),droplet.radius)])
        return sphereNorm

class Vector:
    def __init__(self, x = 0, y = 0, z = 0):
        self.v = [x, y, z]


