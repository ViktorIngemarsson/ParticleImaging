from Vector import Vector
import numpy as np

v = Vector(1, 2, 3, 7, 8, 9)

print("Rotations:")
v.xrotation(55.6).disp()
v.yrotation(55.6).disp()
v.zrotation(55.6).disp()

print("uminus/plus")
v.uminus().disp()
v1 = Vector(4, 5, 6, 10, 11, 12)
v.plus(v1).disp()
v.minus(v1).disp()

print("mtimes/times")
v.mtimes(v1).disp()
print(v.times(v1))

print("rdivide")
v.rdivide(3).disp()

print("versor")
v.versor().disp()
print("topoint")
v.topoint().disp()
print("toline")
v.toline().disp()

print("Plotting")
vPlot = Vector(1,1,1,1,1,1)
vPlot2 = Vector(3,2,5,2,2,4)
vPlot3 = Vector(3,4,4,6,7,1)

vPlot.plot_multiple_vectors([vPlot2,vPlot3])

