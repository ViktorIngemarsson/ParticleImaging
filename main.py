from Point import Point
from Plane import Plane
# def createImage():
# print("This line will be printed.")
#    return 3

# def displayImage():
# print("Should print image")

# a = createImage()

# def refraction(xDim, yDim, dropletRadius, posX, posY, posZ, distZ, rayDensity, refractiveIndexMedium1,


# refraction(2,2,2,0,0,5,10,1,1,1)


# Unittest of plane
p1 = Point(0, 0, 0)
p2 = Point(1, 0, 0)
p3 = Point(0, 1, 0)
pl = Plane(p1, p2, p3)

pl.disp()

pl.xrotation(0.2).disp()

p3 = Point(3, 1, 7)

kl = pl.translate(p3)

kl.disp()
print('jk')
pl.disp()
k = pl.xrotation(0.2)
k.disp()
#l = k.yrotation(0.22)
#o = l.zrotation(0.67)
#o.disp()

