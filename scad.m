function f= scad(x,lam)
f = lam * (x <= lam) + (3.7*lam - x) * (3.7*lam - x > 0) * (x > lam) / 2.7;
end