syms x
syms f(x,y)
syms s
syms t

Rdd = dirac((x.*cos(t) + y.*sin(t)) - s);
Rd  = int(Rdd,'x')
R   = int(Rd,'y')