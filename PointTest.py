from Point import Point
import numpy as np

# TODO: Testa det som innehållet Vector nedan. Rad 17 och 62

#Test-object
p = Point(1, 2, 3)

# Disp
p.disp()

# Translate
dP1 = Point(7, 8, 9)
output = p.translate(dP1)
output.disp()
p.disp()
# = Vector(7, 8, 9, 88, 5, 33)
# output = p.translate(dP2)

# Rotation: x,y,z
print("Rotations:")
p.xrotation(55.6).disp()
p.yrotation(55.6).disp()
p.zrotation(55.6).disp()
# # Numel
print(p.numel())
# Size
print("Size:")
print(p.size())
p2 = Point([1,2],[2,3],[3,4])
print(p2.size())
# uplus, uminus
print("uplus/uminus:")
p.uplus().disp()
p.uminus().disp()
# plus, minus
print("plus/minus:")
p1 = Point(7, 8, 9)
p.plus(p1).disp()
p.minus(p1).disp()
# mtimes, times
print("mtimes/times:")
p.mtimes(p1).disp()
print(p.times(p1))
# rdivide
print("rdivide:")
p.rdivide(4).disp()
# norm
print("Norm:")
print(p.norm())
# normalize
print("normalize:")
p.normalize().disp()
# angle
print("angle:")
print(p.angle(p1))
# # toline
print("toline:")
p.toline()
# # tovector
print("tovector:")
# p.tovector()
print("Plotting")
#pPlot = Point([1,2], [2,3], [3,4])
#p.plot()

