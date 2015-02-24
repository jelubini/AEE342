% Joel Lubinitsky
% AEE 342 - HW6: Incompressible Flow over Airfoils (2)
% 02/25/15

clear all
close all
clc

alpha = 0;
c     = pi;

m  = 2;
p  = 4;
tt = 12;

m  = m * 0.01;
p  = p * 0.1;
tt = tt * 0.01;

zeta     = linspace(pi, 2 * pi, 100);
xCamber  = 0.5 * (1 + cos(zeta));
indexP   = find(xCamber < p, 1, 'last');
xCamber1 = xCamber(1 : indexP);
xCamber2 = xCamber(indexP + 1 : end);
zCamber1 = 0.250 .* (0.800 .* xCamber1 - xCamber1 .^ 2);
zCamber2 = 0.111 .* (0.200 + 0.800 .* xCamber2 - xCamber2 .^ 2);
zCamber  = [zCamber1, zCamber2];

dz1dx = 1 ./ 5 - xCamber1 ./ 2;
dz2dx = 111 ./ 1250 - (111 .* xCamber2) ./ 500;
% zeta = acos(1 - (2 * p) / c)
theta = acos(1 - 2 * p)


% int1 = trapz(

% A0 = 




figure(1)
hold on
axis([-1 2 -1 1])
plot(xCamber, zCamber)