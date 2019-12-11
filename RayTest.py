from Ray import Ray
from Vector import Vector
from Point import Point

v = Vector(1,2,3,4,5,6)
P = 1
pol = Vector(3,1,6,0,1,1)
r = Ray(v,P,pol)

r.disp()

print("Translate")
dP = Point(1,6,8)
dV = Vector(0,9,0,8,7,1)
r.translate(dP).disp()
r.translate(dV).disp()

print("Rotations:")
r.xrotation(55.6).disp()
r.yrotation(55.6).disp()
r.zrotation(55.6).disp()

print("numel")
print(r.numel())

print("plus")
r.uplus().disp()
print("minus")
r.uminus().disp()
print("size")
print(r.size())

print("Angle")
v2 = Vector(7,2,6,3,2,5)
P2 = 1
pol2 = Vector(5,6,1,1,3,3)
r2 = Ray(v2,P2,pol2)
print(r.angle(r2))

print("Versor")
r.versor().disp()

print("ToLine")
r.toline().disp()