from Spherical import Spherical
from Point import Point
from SLine import SLine

p = Point(0, 0, 0)
a = Spherical(p, 2)
a.disp()

dp = Point(1, 0, 0)
b = a.translate(dp)
b.disp()

phi = 45
c = b.xrotation(phi)
c.disp()

d = c.yrotation(phi)
d.disp()

e = d.zrotation(phi)
e.disp()

f = e.numel()
print(f)

v = None
g = e.size(v)
print(g)

#dl = SLine(p, dp)
#nn = 1
#h = e.intersectionpoint(dl, nn)
#h.disp()

p3 = Point(0, 0, 3)
i = e.perpline(p3)