% Joel Lubinitsky
% AEE 342 - HW6 Coefficient Derivations
% 02/24/15


syms alpha z1 z2 theta n x

z1    = 0.250 * (0.800 * x - x ^ 2);
z2    = 0.111 * (0.200 + 0.800 * x - x ^ 2);
dzdx1 = diff(z1, x);
dzdx2 = diff(z2, x);
x     = 0.5 * (1 - cos(theta));
A0    = alpha - (1/pi) * (int(dzdx1, theta, 0, acos(0.2)) + int(dzdx2, theta, acos(0.2), pi))
An    = (2/pi) * (int(dzdx1  * cos(n * theta), theta, 0, acos(0.2)) + int(dzdx2  * cos(n * theta), theta, acos(0.2), pi))
