function f= keru(u,h)
u = u/h;
f=0.75 * (1-u.^2) .* (abs(u) <= 1) /h;
end