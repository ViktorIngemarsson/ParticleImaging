import numpy as np
import cmath

def rayHits(x, y, posX, posY, radius, rayDensity):
    if (np.power(x - posX*rayDensity, 2)+np.power(y-posY*rayDensity, 2)) < np.power(radius*rayDensity, 2):
        return True
    else:
        return False

def getzposonsphere(x,y,radius,posX,posY,posZ):
    relativeZ = np.sqrt(np.power(radius,2) - np.power(posX -x,2) - np.power(posY - y,2))
    return (posZ - relativeZ)

def lengthofvector(v):
    dist = 0
    for i in v:
        dist = dist + np.power(i,2)
    dist = np.sqrt(dist)
    return dist

def snellsLawVector(v1, sphereNorm, refractiveIndexMedium1, refractiveIndexMedium2):
    lengthOfSphereNorm = lengthofvector(sphereNorm)
    lengthOfVector = lengthofvector(v1)

    normalisedsphereNorm = sphereNorm/lengthOfSphereNorm
    normalisedv1 = v1/lengthOfVector

    eta =refractiveIndexMedium1/refractiveIndexMedium2
    sqrtValue = 1-eta*eta*(1-np.power(-normalisedsphereNorm*normalisedv1,2))
    v2 = np.dot(eta,normalisedv1) + (-np.dot(eta, normalisedsphereNorm*normalisedv1) - np.sqrt(1-eta*eta*(1-np.power(-normalisedsphereNorm*normalisedv1,2))))
    return v2

def gerSecondPosOnSphere(x,y,z, vector, posX,posY,posZ,radius):
    a = np.power(vector[0],2) + np.power(vector[1],2) + np.power(vector[2],2)
    b = 2*x*vector[0]-2*vector[0]*posX + 2 * y *vector[1] - 2* vector[1] * posY + 2 * z * vector[2] - 2 * vector[2] * posZ
    c = np.power(x,2) - 2 * x * posX + np.power(posX,2) +np.power(y,2) - 2 * y * posY + np.power(posY,2) + np.power(z,2) - 2 * z * posZ + np.power(posZ,2) - radius*radius
    delta = (b ** 2) - (4 * a * c)
    solution1 = (-b - cmath.sqrt(delta)) / (2 * a)
    solution2 = (-b + cmath.sqrt(delta)) / (2 * a)
    xOut = x + vector[0]*solution2.real
    yOut = y + vector[1]*solution2.real
    zOut = z + vector[2]*solution2.real
    return [xOut,yOut,zOut]

def getSphereNorm(x, y, z, posX, posY, posZ, radius):

    sphereNorm = [(x-posX)/radius,  (y-posY)/radius, (z-posZ)/radius]
    return sphereNorm

def refraction(xDim, yDim, dropletRadius, posX, posY, posZ, distZ, rayDensity, refractiveIndexMedium1, refractiveIndexMedium2):
    image = [0] * xDim * yDim * rayDensity * rayDensity
    for iRayX in range(-xDim,xDim,1):
        for iRayY in range(-yDim,yDim,1):
            if rayHits(iRayX, iRayY, posX, posY, dropletRadius, rayDensity):
                V_1 = [0, 0, 1]
                zPosition = getzposonsphere(iRayX,iRayY,dropletRadius,posX,posY,posZ)
                sphereNorm = getSphereNorm(iRayX,iRayY,zPosition,posX, posY, posZ, dropletRadius)
                V_2 = snellsLawVector(V_1, sphereNorm, refractiveIndexMedium1, refractiveIndexMedium2)
                posOut = gerSecondPosOnSphere(iRayX,iRayY,zPosition,V_2,posX,posY, posZ,dropletRadius)
                sphereNorm2 = getSphereNorm(posOut[0], posOut[1], posOut[2], posX, posY, posZ, dropletRadius)
                for i in range(3):
                    sphereNorm2[i] = -sphereNorm2[i]
                V_3 = snellsLawVector(V_2, sphereNorm2, refractiveIndexMedium2, refractiveIndexMedium1)
                print(V_3)
