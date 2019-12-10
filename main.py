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


# Unittest of plane

p1 = Point(2, 3, 5);
p2 = Point(1, 1, 5);
p3 = Point(2.4, 1, 7);
pl = Plane(p1, p2, p3);

pl.disp()

pl.xrotation(0.2).disp()

p3 = Point(3, 1, 7)

kl = pl.translate(p3)

kl.disp()
print('jk')
pl.disp()
k = pl.xrotation(0.2)
l = k.yrotation(0.22)
o = l.zrotation(0.67)

jil = o.numel()
print(jil)
print(o.size())

p5 = Point(2, 6.2, 2)

L = SLine(p1,p5)

pl.intersectionpoint(L)


